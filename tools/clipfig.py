import os
import matplotlib.pyplot as plt
import pyperclip
import tempfile

def clipfig(fig=None):
    """
    Copies the content of a Matplotlib figure to the clipboard as a vector graphic.

    Parameters:
        fig (matplotlib.figure.Figure): The Matplotlib figure to copy.
    """
    # if no figure is provided, use the current figure
    if fig is None:
        fig = plt.gcf()

    # embed the fonts in the SVG file
    plt.rcParams['svg.fonttype'] = 'none'

    # Save the figure as a temporary SVG file
    with tempfile.NamedTemporaryFile(suffix='.svg', delete=False) as temp_file:
        fig.savefig(temp_file, format='svg')
        temp_file_path = temp_file.name

    # Read the content of the SVG file
    with open(temp_file_path, 'r') as temp_file:
        content = temp_file.read()

    # Copy the content to the clipboard
    pyperclip.copy(content)

    # Delete the temporary SVG file
    os.remove(temp_file_path)

    print('Content copied to clipboard')
