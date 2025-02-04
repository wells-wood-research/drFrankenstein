import numpy as np
import openmm as mm
import parmed
import torch
import torch.nn as nn
import torch.optim as optim


class PositiveFloatConstraint:
    def __call__(self, parameter):
        return torch.clamp(parameter, min=0.0)


class IntegerConstraint:
    def __call__(self, parameter):
        return torch.round(parameter)


class PhaseConstraint:
    def __call__(self, parameter):
        return 180.0 * torch.round(parameter / 180.0)


class StructureModifier(torch.nn.Module):
    def __init__(self, structure: parmed.Structure):
        super().__init__()
        self.structure = structure
        self.phi_k = torch.tensor(
            [dih.phi_k for dih in self.structure.dihedral_types],
        )
        self.per = torch.tensor([dih.per for dih in self.structure.dihedral_types])

        self.phase = torch.tensor([dih.phase for dih in self.structure.dihedral_types])

        self.phi_k_phase = nn.Parameter(
            torch.mul(
                self.phi_k,
                torch.where(self.phase > 179, -1, 1),
            )
        )
        self.K = nn.Parameter(torch.tensor([0.00]))
        self.K.requires_grad = True
        self.phi_k_phase.requires_grad = True

    def deriv(self, xyz):
        """
        Here aproximated this by taking the derivative of the truncated Fourrrier series
        """

        return torch.mul(
            torch.sin(self.per), torch.sin(torch.mul(self.per, self.phi_k_phase))
        )

    def forward(
        self,
        xyz,
    ):
        self.update_structure()

        # current_MM_energies = self.compute_MM_energy(xyz)
        deriv = self.deriv(xyz)

        current_MM_energies = self.compute_MM_energy(xyz)

        return Custom_Function.apply(
            current_MM_energies,
            deriv,
            self.phi_k_phase,
            self.K,
        )
        # Update dihedral parameters

    def update_structure(self):
        self.phi_k = torch.abs(self.phi_k_phase)
        self.phase = torch.where(self.phi_k_phase < 0, 180.0, 0.0)

        for dih, phi_k, phase, per in zip(
            self.structure.dihedral_types, self.phi_k, self.phase, self.per
        ):
            dih.phi_k = float(phi_k)
            dih.per = int(abs(per))
            dih.phase = float(phase)
        self.write_frcmod()

    def compute_MM_energy(self, frames):
        parm = self.structure
        system = parm.createSystem(nonbondedMethod=None)

        # Make the context and set the positions
        context = mm.Context(system, mm.VerletIntegrator(0.001))

        energy_total = []

        for positions in frames:
            context.setPositions(positions)
            # Find the energy decomposition
            energy_total.append(
                parmed.openmm.energy_decomposition(parm, context)["total"] + self.K[0]
            )

        return torch.tensor(energy_total, requires_grad=True)

    def write_frcmod(self, file_out="test.frcmod"):
        parmed.tools.writeFrcmod(self.structure, file_out).execute()


class Custom_Function(torch.autograd.Function):
    @staticmethod
    def forward(
        ctx,
        MM_energies,
        deriv,
        phi_k_phase,
        K,
    ):
        # Save the inputs for the backward pass
        print("forward")
        ctx.save_for_backward(K, phi_k_phase, deriv)
        return MM_energies + K

    @staticmethod
    def backward(ctx, grad_output):
        """
        In the backward pass we receive a Tensor containing the gradient of the loss
        with respect to the output, and we need to compute the gradient of the loss
        with respect to the input.
        Here aproximated this by taking the derivative of the truncated Fourrrier series
        """
        print("entering backward")
        # Retrieve the saved inputs
        (K, phi_k_phase, deriv) = ctx.saved_tensors
        print(grad_output)
        grad_K = torch.tensor([grad_output.mean() - K])
        print(grad_K)
        if grad_K > 50:
            grad_phi_k_phase = torch.zeros(len(deriv.t()))
        else:
            grad_phi_k_phase = grad_K * 0.005 * deriv.t() * phi_k_phase.t()

        return (
            None,
            None,
            grad_phi_k_phase,
            grad_K,
        )


def train_model(structure):
    structure_modifier = StructureModifier(structure)
    optimizer = torch.optim.Adam(structure_modifier.parameters(), lr=learning_rate)

    for epoch in range(num_epochs):
        for i, (positions, target_MM_energies) in enumerate(train_loader):
            optimizer.zero_grad()
            predicted_MM_energies = structure_modifier.compute_MM_energy(positions)
            loss = F.mse_loss(predicted_MM_energies, target_MM_energies)
            loss.backward()
            optimizer.step()

    import torch.optim as optim

    # Define your input and target data
    input_data = ...  # list of energy values that you want to fit to
    target_data = ...  # list of corresponding atom positions

    # Define your model and optimizer
    model = StructureModifier(structure)
    optimizer = optim.Adam(model.parameters(), lr=0.001)

    # Define your loss function
    loss_fn = nn.MSELoss()

    # Train the model
    for epoch in range(num_epochs):
        running_loss = 0.0
        for i in range(len(input_data)):
            # Zero the parameter gradients
            optimizer.zero_grad()

            # Forward pass
            output = model(target_data[i])

            # Compute the loss
            loss = loss_fn(output, input_data[i])

            # Backward pass and optimize
            loss.backward()
            print(loss)
            optimizer.step()
            print(model.parameters())
            # Print statistics
            running_loss += loss.item()
            if i % 100 == 99:
                print(f"Epoch {epoch+1}, batch {i+1}: loss={running_loss/100:.4f}")
                running_loss = 0.0

    structure_modifier.writefrcmod(file_out="optimised dihedrals.frcmod").execute()


import time


def train(model, optimizer, loss_fn, train_dl, val_dl, epochs=100, device="cpu"):
    print(
        "train() called: model=%s, opt=%s(lr=%f), epochs=%d, device=%s\n"
        % (
            type(model).__name__,
            type(optimizer).__name__,
            optimizer.param_groups[0]["lr"],
            epochs,
            device,
        )
    )

    history = {}  # Collects per-epoch loss and acc like Keras' fit().
    history["loss"] = []
    history["val_loss"] = []
    history["acc"] = []
    history["val_acc"] = []

    start_time_sec = time.time()
    import sys

    for epoch in range(1, epochs + 1):
        # --- TRAIN AND EVALUATE ON TRAINING SET -----------------------------
        # Set the model to training mode
        model.train()
        print("heho?? ")
        train_loss = 0.0
        num_train_correct = 0
        num_train_examples = 0

        for batch in train_dl:
            optimizer.zero_grad()

            x = batch[0].numpy()
            y = batch[1].to(device)
            print("Starting")

            yhat = model(x)
            print(yhat, y)

            loss = loss_fn(yhat, y)
            print(f"Loss : {loss}")

            loss.backward()
            for param in model.parameters():
                print(f" gradient: {param.grad}")
                if param.grad is None:
                    sys.exit()
            print(list(model.parameters()))
            optimizer.step()
            print(list(model.parameters()))
            train_loss += loss.data.item() * y.size(0)

        # --- EVALUATE ON VALIDATION SET -------------------------------------
        model.eval()
        val_loss = 0.0
        num_val_correct = 0
        num_val_examples = 0

        for batch in val_dl:
            x = batch[0].numpy()
            y = batch[1].to(device)
            yhat = model(x)
            loss = loss_fn(yhat, y)

            val_loss += loss.data.item() * y.size(0)

        if epoch == 1 or epoch % 10 == 0:
            print(
                "Epoch %3d/%3d, train loss: %5.2f,  val loss: %5.2f, "
                % (
                    epoch,
                    epochs,
                    train_loss,
                    val_loss,
                )
            )

        history["loss"].append(train_loss)
        history["val_loss"].append(val_loss)

    # END OF TRAINING LOOP

    end_time_sec = time.time()
    total_time_sec = end_time_sec - start_time_sec
    time_per_epoch_sec = total_time_sec / epochs
    print()
    print("Time total:     %5.2f sec" % (total_time_sec))
    print("Time per epoch: %5.2f sec" % (time_per_epoch_sec))
    model.write_frcmod("final.frcmod")
    return history


if __name__ == "__main__":
    import pickle as pkl

    data = pkl.load(open("5_points.pkl", "rb"))

    x_coordinates = [d[1] for d in data]
    energies = [d[0] for d in data]

    train_dl_x = torch.tensor(x_coordinates)
    train_dl_y = torch.tensor(energies)
    train_dl = [(train_dl_x, train_dl_y)]
    structure = parmed.load_file("mol.prmtop")

    model = StructureModifier(structure)

    optimizer = torch.optim.SGD(
        model.parameters(), lr=0.1
    )  # optim.Adam(model.parameters(), lr=0.01)
    loss_fn = nn.MSELoss()
    loss_fn = torch.nn.functional.mse_loss

    model.write_frcmod()
    train(model, optimizer, loss_fn, train_dl, train_dl, epochs=5, device="cpu")
