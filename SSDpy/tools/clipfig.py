import os
import matplotlib.pyplot as plt
import pyperclip
import tempfile
import json

import subprocess

def clipfig(fig=None, save=False):
    """
    Copies the content of a Matplotlib figure to the clipboard as a vector graphic.

    Parameters:
        fig (matplotlib.figure.Figure): The Matplotlib figure to copy.
    """
    # if no figure is provided, use the current figure
    if fig is None:
        fig = plt.gcf()

    if not save:
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
    
    else:
        # Open Finder at the folder containing the file
        folder_path = "/tmp/ssd"
        os.makedirs(folder_path, exist_ok=True)

        # Optionally save a file in it
        file_path = os.path.join(folder_path, "clipfig_00001.pdf")
        # check if the file already exists
        i = 1
        while os.path.exists(file_path):
            file_path = os.path.join(folder_path, f"clipfig_{i:05d}.pdf")
            i += 1
            if i > 99999:
                raise ValueError("Too many files, please clean up the folder")
        # Save the figure as a PDF file
        fig.savefig(file_path, format="pdf", bbox_inches="tight", pad_inches=0)

        # Open Finder and select the file
        subprocess.run(["open", "-R", file_path])

        print("Temporary file saved at:", file_path)

import matplotlib.pyplot as plt

def extract_figure_data(fig=None):
    """
    Extracts all artist information from a matplotlib figure and
    stores it in a structured dictionary.
    """
    if fig is None:
        return

    fig_data = {"axes": []}
    
    for ax in fig.get_axes():
        # NEED TO BE UPDATED
        ax_info = {
            "title": ax.get_title(),
            "xlabel": ax.get_xlabel(),
            "ylabel": ax.get_ylabel(),
            "lines": [],
            "patches": [],
            "collections": [],
            "texts": [],
            "images": []
        }
        
        # Lines (Line2D)
        for line in ax.get_lines():
            ax_info["lines"].append({
                "x": line.get_xdata().tolist(),
                "y": line.get_ydata().tolist(),
                "color": line.get_color(),
                "linestyle": line.get_linestyle(),
                "linewidth": line.get_linewidth(),
                "marker": line.get_marker(),
                "label": line.get_label()
            })
        
        # Patches (Rectangle, Circle, etc.)
        for patch in ax.patches:
            ax_info["patches"].append({
                "type": type(patch).__name__,
                "xy": patch.get_xy() if hasattr(patch, "get_xy") else None,
                "width": patch.get_width() if hasattr(patch, "get_width") else None,
                "height": patch.get_height() if hasattr(patch, "get_height") else None,
                "radius": patch.get_radius() if hasattr(patch, "get_radius") else None,
                "facecolor": patch.get_facecolor(),
                "edgecolor": patch.get_edgecolor()
            })
        
        # Collections (e.g., scatter, pcolormesh, contourf)
        for coll in ax.collections:
            entry = {"type": type(coll).__name__}
            if hasattr(coll, "get_offsets"):  # e.g. scatter
                entry["offsets"] = coll.get_offsets().tolist()
            if hasattr(coll, "get_array"):    # e.g. pcolormesh, contourf
                arr = coll.get_array()
                entry["array"] = arr.tolist() if arr is not None else None
            if hasattr(coll, "get_facecolors"):
                entry["facecolors"] = coll.get_facecolors().tolist()
            ax_info["collections"].append(entry)
        
        # Text (annotations, labels)
        for text in ax.texts:
            ax_info["texts"].append({
                "string": text.get_text(),
                "position": text.get_position(),
                "fontsize": text.get_fontsize(),
                "color": text.get_color()
            })
        
        # Images (imshow)
        for img in ax.images:
            ax_info["images"].append({
                "array": img.get_array().tolist(),
                "extent": img.get_extent(),
                "cmap": img.get_cmap().name,
                "clim": img.get_clim()
            })
        
        fig_data["axes"].append(ax_info)
    
        # Save figure data in a temporary file
                # Open Finder at the folder containing the file
        folder_path = "/tmp/ssd"
        os.makedirs(folder_path, exist_ok=True)

        # Optionally save a file in it
        file_path = os.path.join(folder_path, "fig_data_00001.pdf")
        # check if the file already exists
        i = 1
        while os.path.exists(file_path):
            file_path = os.path.join(folder_path, f"fig_data_{i:05d}.pdf")
            i += 1
            if i > 99999:
                raise ValueError("Too many files, please clean up the folder")
        # Save the dictionary
        with open(file_path, "w") as f:
            json.dump(fig_data, f)
        print("Figure data saved at:", file_path)
    return file_path, fig_data

import matplotlib.pyplot as plt
import numpy as np

def recreate_figure(fig_data_file, fig=None):
    """
    Rebuilds a matplotlib figure from the extracted figure data dictionary.
    """
    folder_path = "/tmp/ssd"
    os.makedirs(folder_path, exist_ok=True)
    # Check if the file already exists
    file_path = os.path.join(folder_path, fig_data_file)
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} does not exist.")

    with open(file_path, "r") as f:
        fig_data = json.load(f)
        
    # Create a new figure
    if fig is None:
        fig, axes_list = plt.subplots(
            nrows=len(fig_data["axes"]),
            ncols=1,
            squeeze=False
        )
        axes_list = axes_list.flatten()
    else:
        axes_list = fig.get_axes()
        if len(axes_list) != len(fig_data["axes"]):
            raise ValueError("The number of axes in the figure does not match the data.")
        
    for ax, ax_info in zip(axes_list, fig_data["axes"]):
        # Titles and labels
        ax.set_title(ax_info["title"])
        ax.set_xlabel(ax_info["xlabel"])
        ax.set_ylabel(ax_info["ylabel"])
        
        # Lines
        for line_info in ax_info["lines"]:
            ax.plot(
                line_info["x"], 
                line_info["y"], 
                color=line_info["color"],
                linestyle=line_info["linestyle"],
                linewidth=line_info["linewidth"],
                marker=line_info["marker"],
                label=line_info["label"] if line_info["label"] != "_nolegend_" else None
            )
        
        # Patches
        for patch_info in ax_info["patches"]:
            t = patch_info["type"]
            if t == "Rectangle":
                from matplotlib.patches import Rectangle
                ax.add_patch(Rectangle(
                    patch_info["xy"],
                    patch_info["width"],
                    patch_info["height"],
                    facecolor=patch_info["facecolor"],
                    edgecolor=patch_info["edgecolor"]
                ))
            elif t == "Circle":
                from matplotlib.patches import Circle
                ax.add_patch(Circle(
                    patch_info["xy"],
                    patch_info["radius"],
                    facecolor=patch_info["facecolor"],
                    edgecolor=patch_info["edgecolor"]
                ))
            # Add more patch types as needed
        
        # Collections (scatter, pcolormesh, contourf)
        for coll_info in ax_info["collections"]:
            t = coll_info["type"]
            if t == "PathCollection" and "offsets" in coll_info:
                offsets = np.array(coll_info["offsets"])
                facecolors = np.array(coll_info.get("facecolors", []))
                # If you stored scalar array for coloring
                if "array" in coll_info and coll_info["array"] is not None:
                    c = np.array(coll_info["array"])
                else:
                    c = facecolors
                ax.scatter(offsets[:,0], offsets[:,1], c=c)
            elif t == "QuadMesh" and "array" in coll_info:
                # Here you'd need the original x/y grid if you want exact reconstruction
                z = np.array(coll_info["array"])
                ax.pcolormesh(z.reshape(int(np.sqrt(z.size)), -1))
            # Other collection types (PolyCollection, etc.) can be added similarly
        
        # Texts
        for text_info in ax_info["texts"]:
            ax.text(
                *text_info["position"],
                text_info["string"],
                fontsize=text_info["fontsize"],
                color=text_info["color"]
            )
        
        # Images
        for img_info in ax_info["images"]:
            arr = np.array(img_info["array"])
            ax.imshow(
                arr, 
                extent=img_info["extent"],
                cmap=plt.get_cmap(img_info["cmap"]),
                clim=img_info["clim"],
                origin='upper'
            )
    
    fig.tight_layout()
    return fig

def duplfig(fig1, fig2):
    """
    Duplicates the content of one Matplotlib figure to another.
    
    Parameters:
        fig1 (matplotlib.figure.Figure): The source figure to copy from.
        fig2 (matplotlib.figure.Figure): The target figure to copy to. If None, a new figure is created.
    """
    if fig2 is None or not isinstance(fig2, plt.Figure):
        fig2 = plt.figure()
    
    file_path, fig_data = extract_figure_data(plt.figure(fig1))
    recreate_figure(file_path, fig=fig2)
    
    for ax1 in fig1.get_axes():
        ax2 = fig2.add_subplot(111)  # Add a new subplot
        ax2.set_title(ax1.get_title())
        ax2.set_xlabel(ax1.get_xlabel())
        ax2.set_ylabel(ax1.get_ylabel())
    
    return fig2