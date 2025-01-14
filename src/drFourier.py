import os
from os import path as p
import pandas as pd
import numpy as np
import plotly.graph_objects as go
# Import the required packages
from scipy.fft import fft, rfft
from scipy.fft import fftfreq, rfftfreq
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt

def dummy_inputs():
    dataDir = "/home/esp/scriptDevelopment/drFrankenstein/ALA_outputs/torsion_C-CA-N-CC1/scan_data"
    energyCsv = p.join(dataDir, "final_scan_energies.csv")
    energyDf = pd.read_csv(energyCsv)

    torsionTag = 'C-CA-N-CC1'
    return dataDir, energyDf, torsionTag

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲

def main():
    
    dataDir, energyDf, torsionTag = dummy_inputs()

    energyData = energyDf[torsionTag].values
    angleData = energyDf['Angle'].values


    energyData_padded = np.concatenate((energyData, energyData[::-1], energyData, energyData[::-1], energyData, energyData[::-1]))
    angleData_padded = np.concatenate((angleData, angleData + 360, angleData + 720, angleData + 1080, angleData + 1440, angleData + 1800))

    fourier = Fourier(energyData_padded, angleAxis=  angleData_padded, sampling_rate=0.1, out_file="fourier.html")


    fourier.plot_angle_frequency(t_title="ECG Signal", f_title="ECG Spectrum",
                                 t_ylabel="Amplitude[mV]")
    
    fourierData = fourier.extract_cosine_parameters()
    print(fourierData)


    cosineDf = fourier.cosine_functions()
    print(cosineDf)


    reconstructedEnergy = fourier.reconstruct_signal(3)
    plot_interactive_lines(angleData_padded, cosineDf, energyData_padded, reconstructedEnergy, nComponants=4, fileName = "plot.html")






def create_amber_torsion_params(torsionTag, fourierData, nComponants):
    pass

#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲


# Building a class Fourier for better use of Fourier Analysis.
class Fourier:
    """
    Apply the Discrete Fourier Transform (DFT) on the signal using the Fast Fourier 
    Transform (FFT) from the scipy package.

    Example:
        fourier = Fourier(signal, sampling_rate=2000.0, out_file='output.html')
    """

    def __init__(self, signal, angleAxis, sampling_rate, out_file):
        """
        Initialize the Fourier class.

        Args:
            signal (np.ndarray): The samples of the signal
            sampling_rate (float): The sampling per second of the signal
            out_file (str): The file name to save the output graph
        """
        self.signal = signal
        self.sampling_rate = sampling_rate
        self.out_file = out_file
        self.angle_step = 1.0 / self.sampling_rate
        self.duration = len(self.signal) / self.sampling_rate
        self.angle_axis = angleAxis
        self.frequencies = rfftfreq(len(self.signal), d=self.angle_step)
        self.fourier = rfft(self.signal)

    def amplitude(self):
        """
        Method of Fourier

        Returns:
            numpy.ndarray of the actual amplitudes of the sinusoids.
        """
        return 2 * np.abs(self.fourier) / len(self.signal)

    def phase(self, degree=True):
        """
        Method of Fourier

        Args:
            degree: To choose the type of phase representation (Radian, Degree).
                    By default, it's in radian. 

        Returns:
            numpy.ndarray of the phase information of the Fourier output.
        """
        return np.angle(self.fourier, deg=degree)

    def plot_spectrum(self, interactive=False):
        """
        Plot the Spectrum (Frequency Domain) of the signal either using the matplotlib
        package, or plot it interactive using the plotly package.

        Args:
            interactive: To choose if you want the plot interactive (True), or not
            (False). The default is the spectrum non-interactive.

        Returns:
            A plot of the spectrum.
        """
        if interactive:
            trace = go.Scatter(x=self.frequencies, y=self.amplitude())
            data = [trace]
            layout = go.Layout(title=dict(text='Spectrum',
                                          x=0.5,
                                          xanchor='center',
                                          yanchor='top',
                                          font=dict(size=25, family='Arial, bold')),
                               xaxis=dict(title='Frequency[Hz]'),
                               yaxis=dict(title='Amplitude'))
            fig = go.Figure(data=data, layout=layout)
            fig.write_html(self.out_file)
        else:
            plt.figure(figsize=(10, 6))
            plt.plot(self.frequencies, self.amplitude())
            plt.title('Spectrum')
            plt.ylabel('Amplitude')
            plt.xlabel('Frequency[Hz]')
            plt.savefig(self.out_file)

    def plot_angle_frequency(self, t_ylabel="Amplitude", f_ylabel="Amplitude",
                            t_title="Signal (angle Domain)",
                            f_title="Spectrum (Frequency Domain)"):
        """
        Plot the Signal in angle Domain and Frequency Domain using plotly.

        Args:
            t_ylabel (String): Label of the y-axis in angle-Domain
            f_ylabel (String): Label of the y-axis in Frequency-Domain
            t_title (String): Title of the angle-Domain plot
            f_title (String): Title of the Frequency-Domain plot 

        Returns:
            Two figures: the first is the angle-domain, and the second is the
                         frequency-domain.
        """
        # The Signal (Angle-Domain)
        angle_trace = go.Scatter(x=self.angle_axis, y=self.signal)
        angle_domain = [angle_trace]
        layout = go.Layout(title=dict(text=t_title,
                                      x=0.5,
                                      xanchor='center',
                                      yanchor='top',
                                      font=dict(size=25, family='Arial, bold')),
                           xaxis=dict(title='angle[degrees]'),
                           yaxis=dict(title=t_ylabel),
                           width=1000,
                           height=400)
        fig = go.Figure(data=angle_domain, layout=layout)
        fig.write_html(self.out_file.replace('.html', '_angle.html'))

        # The Spectrum (Frequency-Domain)
        freq_trace = go.Scatter(x=self.frequencies, y=self.amplitude())
        frequency_domain = [freq_trace]
        layout = go.Layout(title=dict(text=f_title,
                                      x=0.5,
                                      xanchor='center',
                                      yanchor='top',
                                      font=dict(size=25, family='Arial, bold')),
                           xaxis=dict(title='Frequency[Hz]'),
                           yaxis=dict(title=f_ylabel),
                           width=1000,
                           height=400)
        fig = go.Figure(data=frequency_domain, layout=layout)
        fig.write_html(self.out_file.replace('.html', '_freq.html'))
    def extract_cosine_parameters(self):
        """
        Extracts the amplitudes, frequencies, and phases for each cosine function
        from the Fourier Transform, sorted by amplitude in descending order.

        Returns:
            pd.DataFrame: A DataFrame containing the amplitudes, frequencies,
                          and phases of the Fourier components, sorted by amplitude.
        """
        amplitudes = self.amplitude()
        phases = self.phase()
        
        # Get the indices that would sort the amplitudes in descending order
        sorted_indices = np.argsort(-amplitudes)
        
        # Sort the data by these indices
        data = {
            'Frequency ': self.frequencies[sorted_indices],
            'Amplitude': amplitudes[sorted_indices],
            'Phase (Degrees)': phases[sorted_indices]
        }
        
        return pd.DataFrame(data)
    
    def cosine_functions(self):
        """
        Constructs cosine functions from the Fourier Transform components,
        ordered by amplitude, and includes the Angle column.

        Returns:
            pd.DataFrame: A DataFrame where each column represents a cosine
                          function corresponding to a Fourier component,
                          ordered by amplitude, with an Angle column.
        """
        amplitudes = self.amplitude()
        phases = self.phase()

        # Sort indices by amplitude in descending order
        sorted_indices = np.argsort(-amplitudes)

        cosine_data = {}

        # Construct each cosine function in order of amplitude
        for i in sorted_indices:
            amplitude = amplitudes[i]
            frequency = self.frequencies[i]
            phase = phases[i]
            cosine_data[f'cosine_{i}'] = amplitude * np.cos(
                2 * np.pi * frequency * self.angle_axis + phase
            )

        # Create DataFrame
        cosine_df = pd.DataFrame(cosine_data)
        cosine_df['Angle'] = self.angle_axis  # Assuming Angle corresponds to angle_axis

        return cosine_df
    
    def reconstruct_signal(self, n_components):
        """
        Reconstructs the energy signal using the first N cosine functions
        ordered by amplitude.

        Args:
            n_components (int): The number of cosine functions to use for
                                reconstruction.

        Returns:
            np.ndarray: The reconstructed signal.
        """
        amplitudes = self.amplitude()
        phases = self.phase()

        # Sort indices by amplitude in descending order
        sorted_indices = np.argsort(-amplitudes)

        # Initialize the reconstructed signal
        reconstructed_signal = np.zeros_like(self.angle_axis)

        # Sum the first N cosine functions
        for i in sorted_indices[:n_components]:
            amplitude = amplitudes[i]
            frequency = self.frequencies[i]
            # if frequency == 0:
            #     continue
            phase = phases[i]
            reconstructed_signal += amplitude * np.cos(
                2 * np.pi * frequency * self.angle_axis + phase
            )

        return reconstructed_signal

        
#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲
def plot_interactive_lines(angleData, cosineDf, energyData, reconstructedEnergy, nComponants, fileName='plot.html'):
    # Create a figure
    fig = go.Figure()

    # Add cosine functions as dotted lines
    for col in cosineDf.columns[:-1][:nComponants]:  # Exclude 'Angle' column
        fig.add_trace(go.Scatter(
            x=angleData,
            y=cosineDf[col],
            mode='lines',
            line=dict(dash='dot'),
            name=col
        ))

    # Add energy data as a solid line
    fig.add_trace(go.Scatter(
        x=angleData,
        y=energyData,
        mode='lines',
        line=dict(dash='solid'),
        name='Energy'
    ))

    # Add reconstructed energies as a solid line
    fig.add_trace(go.Scatter(
        x=angleData,
        y=reconstructedEnergy,
        mode='lines',
        line=dict(dash='solid'),
        name='Reconstructed Energy'
    ))

    # Update layout for interactivity
    fig.update_layout(
        title='Interactive Line Plot',
        xaxis_title='Angle',
        yaxis_title='Value',
        legend_title='Components',
        hovermode='x unified'
    )

    # Save the plot as an HTML file
    fig.write_html(fileName)



    

def decompose_energies(energy_df):
    angles = energy_df['Angle'].values
    energies = energy_df['C-CA-N-CC1'].values


    anglesExtended = np.concatenate((angles, angles[::-1], angles, angles[::-1]))
    energiesExtended = np.concatenate((energies, energies[::-1], energies, energies[::-1]))


    print(len(anglesExtended))
    print(anglesExtended)
    print(len(energiesExtended))

    # Perform Fourier Transform
    fourierCoefficients = np.fft.fft(energiesExtended)

    # Get frequencies
    frequencies = np.fft.fftfreq(len(energiesExtended), d=(anglesExtended[1] - anglesExtended[0]))

    # Extract amplitudes and phases
    amplitudes = np.abs(fourierCoefficients) / len(energiesExtended)
    phases = np.angle(fourierCoefficients)

    # Prepare list of cosine parameters
    cosineParams = []
    cosineFunctions = np.zeros((len(energiesExtended), len(energiesExtended)))

    for i, (amplitude, phase) in enumerate(zip(amplitudes, phases)):
        cosineParams.append({
            'amplitude': amplitude,
            'frequency': frequencies[i],
            'phase': phase
        })
        cosineFunctions[:, i] = amplitude * np.cos(
            frequencies[i] * anglesExtended + phase
        )

    # Sort cosine parameters by amplitude in descending order
    sortedIndices = np.argsort(-amplitudes)
    cosineParams = [cosineParams[i] for i in sortedIndices]
    cosineFunctions = cosineFunctions[:, sortedIndices]

    # Create DataFrame for cosine functions
    cosineDf = pd.DataFrame(
        cosineFunctions, 
        columns=[f'cosine_{i}' for i in sortedIndices]
    )
    cosineDf['Angle'] = anglesExtended


    reconstructedEnergies = np.sum(cosineFunctions, axis=1)

    return cosineParams, cosineDf, reconstructedEnergies


#🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲🗲



if __name__ == "__main__":
    main()