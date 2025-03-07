import os
import matplotlib.pyplot as plt


def saveas(folder_path, filename, fig=None):
    """
    Saves a Matplotlib figure as a PDF file with vector graphics.
    This is a function of the SSJPy package.

    Parameters:
        fig (matplotlib.figure.Figure): The Matplotlib figure to save.
        folder_path (str): The path to the directory where the PDF file should be saved.
        filename (str): The name of the PDF file.
    """

    # If no figure is provided, use the current figure
    if fig is None:
        fig = plt.gcf()
        
    # Ensure that the directory exists
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    
    # Set pdf.fonttype to 'embedded'
    plt.rcParams['pdf.fonttype'] = 42

    # Save the figure as a PDF file
    file_path = os.path.join(folder_path, filename)
    fig.savefig(file_path, format='pdf')
