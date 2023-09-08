
import EasySQTools as tools
import matplotlib.pyplot as plt


# displays a color file as a palette
def displayPalette(color_file, markerSize=10):
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


# takes a list and multiplies the number of entries by the size
# l1 = ["item1", "item2", "item3"]
# size = 2
# l2 = ['item1', 'item1', 'item2', 'item2', 'item3', 'item3']
def boxifyList(l1, size):
    l2 = []
    for item in l1:
        for i in range(size):
            l2.append(item)

    return l2


if __name__ == "__main__":
    l1 = ["item1", "item2", "item3"]
    print(boxifyList(l1, 2))

    file = "colors/leiden_color_set_1_gradient.csv"
    # file = "colors/leiden_color_set_1_random.csv"
    # file = "colors/leiden_color_set_2_random.csv"
    # file = "colors/leiden_color_set_3_random.csv"
    # file = "colors/leiden_generated_random_1.txt"
    # file = "colors/leiden_generated_random_2.txt"
    # file = "colors/leiden_generated_random_3.txt"
    displayPalette(file)

