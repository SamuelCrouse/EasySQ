# Sam Crouse
# scrouse2@uwyo.edu

import EasySQTools as tools
import matplotlib.pyplot as plt


def displayPalette(color_file, markerSize=10):
    """
        Documentation

        displayPalette() docs
        =============================
        ----

        * display a list of hex colors in a plot to visualize the palette

        ----

        Parameters:
        ----------
        * color_file="path/to/color_file.txt": Path to colors starting at current working directory.
        * markerSize=10: Size of the markers to be displayed.

        Returns:
        -------
        * None: Shows the plot and returns nothing.

        Examples:
        --------
        1. Calling displayPalette():
            * displayPalette("colors/leiden_color_set_1_random.csv", markerSize=5)
    """

    print("got here 3: {}".format(color_file))
    colors = tools.getColors(color_file=color_file)
    print("colors: {}".format(colors))

    paletteLength = len(colors)
    newColors = boxifyList(colors, size=paletteLength)

    # create the x locations
    # x stays the same till we get to a new color
    x = []
    step = 0
    for i in range(len(newColors)):
        if i % paletteLength == 0:
            step += 1

        x.append(step)

    # create the y locations
    # y increases until we get to a new color, then it resets
    y = []
    for i in range(int(len(newColors) / paletteLength)):
        for j in range(paletteLength):
            y.append(j)

    plt.scatter(x=x, y=y, c=newColors, s=markerSize)

    plt.show()

    return None


def boxifyList(l1, size):
    """
        Documentation

        boxifyList() docs
        =============================
        ----

        * Takes in a list (l1) and multiplies the number of entries by the size (size). Returns l2.

        ----

        Parameters:
        ----------
        * l1=["item1", "item2"]: List of items of any type.
        * size=10: The number of times to multiply and add each item.

        Returns:
        -------
        * list: l2: List of input items multiplied by size.

        Examples:
        --------
        1. Calling boxifyList()
            * l1 = ["item1", "item2", "item3"]
            * size = 2
            * result = boxifyList(l1, size)
            * result = l2 = ['item1', 'item1', 'item2', 'item2', 'item3', 'item3']
    """

    l2 = []
    for item in l1:
        for i in range(size):
            l2.append(item)

    return l2


# testing for DisplayPalette.py
if __name__ == "__main__":
    l1 = ["item1", "item2", "item3"]
    print(boxifyList(l1, 2), "\n")

    # tools.createPalette(200, True, True)

    file = "colors/leiden_color_set_2_gradient.csv"
    # file = "colors/leiden_color_set_1_random.csv"
    # file = "colors/leiden_color_set_2_random.csv"
    # file = "colors/leiden_color_set_3_random.csv"
    # file = "colors/leiden_generated_random_1.txt"
    # file = "colors/leiden_generated_random_2.txt"
    # file = "colors/leiden_generated_random_7.txt"

    displayPalette(file, markerSize=5)
