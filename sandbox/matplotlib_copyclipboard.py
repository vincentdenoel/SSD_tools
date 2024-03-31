from SSDpy.tools.tools import save_figure_as_pdf as ssdsave
import matplotlib.pyplot as plt

# Generate a Matplotlib figure
fig, ax = plt.subplots()
ax.plot([1, 2, 3, 4], [1, 4, 2, 3])

# Save the figure as a PDF file in the "output" folder
ssdsave(fig, '/Users/vincentDenoel/Dropbox/009_Maison/Epslog/_HALLIBURTON/MATLAB CODES/drill_the_stand/','my_figure.pdf')
