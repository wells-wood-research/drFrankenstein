
import numpy as np
import openmm as mm
import parmed
import torch
import torch.nn as nn
import torch.optim as optim


class FourierSeries(nn.Module):
    def __init__(self, order_list):
        super(FourierSeries, self).__init__()
        self.order_list = order_list

        self.a: list = nn.Parameter(torch.randn(len(self.order_list)))
        # self.phi :bool  = nn.Parameter(torch.randn(1))

    def forward(self, dih_angle):
        series = 0

        for i, term_order in enumerate(self.order_list):
            # alternate between 0 and pi for phase
            phi_value = 0 if i % 2 == 0 else np.pi

            term = abs(self.a[i]) * (
                1 + torch.cos(term_order * 2 * np.pi * dih_angle + phi_value)
            )
            series += term

        return series


class SumFourierSeries(nn.Module):
    def __init__(self, num_dihedral, order_lists):
        super(SumFourierSeries, self).__init__()
        self.num_dihedral = num_dihedral
        self.fourier_series = nn.ModuleList(
            [FourierSeries(order_list) for order_list in order_lists]
        )
        self.linear = nn.Linear(self.num_inputs, 1)
        self.K = nn.Parameter(torch.randn(1))

    def forward(self, dih_angles_list):
        series_list = []
        for i in range(self.num_dihedral):
            series_list.append(self.fourier_series[i](dih_angles_list[i]))
        series_sum = self.K + torch.stack(series_list, dim=1).sum(dim=1)
        output = self.linear(series_sum)
        return output


def optimiser():
    input_size = [100, 200, 300]  # length of each input sequence
    max_num_terms = 30  # maximum number of terms to use
    model = SumFourierSeries(input_size, max_num_terms)
    optimizer = optim.Adam(model.parameters(), lr=0.01)
    criterion = nn.MSELoss()
    num_epochs = 100
    # training loop
    for epoch in range(num_epochs):
        optimizer.zero_grad()
        outputs = model(inputs_list)
        loss = criterion(outputs, targets)
        loss.backward()
        optimizer.step()

        # print the loss every 10 epochs
        if epoch % 10 == 0:
            print("Epoch {}, Loss: {:.4f}".format(epoch, loss.item()))
