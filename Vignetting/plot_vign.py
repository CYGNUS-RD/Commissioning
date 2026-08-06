## py .\plot_vig.py .\VIG7.DAT

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import ROOT

def process_file(file_path):
    try:
        with open(file_path, 'r') as file:
            content = file.read()

            # Split the content by tabs and newline characters
            numbers = [value for group in content.split('\n') for value in group.split('\t')]

            # Remove empty strings from the list
            numbers = list(filter(lambda x: x != '', numbers))

            count = len(numbers)

            if count > 0:
                # Convert the numbers to float
                numbers = list(map(float, numbers))

                # Define the grid dimensions
                width = 4096
                height = 2304

                # Check if the number of elements matches the grid size
                if count != width * height:
                    print(f"Warning: The number of elements ({count}) does not match the grid size ({width}x{height}={width*height}).")
                    # Pad with zeros if there are fewer numbers, or truncate if there are too many
                    if count < width * height:
                        numbers.extend([0.0] * (width * height - count))
                    else:
                        numbers = numbers[:width * height]

                # Count how many numbers are equal to 0
                zero_count = numbers.count(0.0)
                print(f"Number of zeros in the data: {zero_count}")

                # Reshape the list into the specified grid size
                data = np.array(numbers).reshape((height, width))

                #save in root file
                stem, _ = os.path.splitext(file_path)
                outname = stem+'.root' 
                rf = ROOT.TFile(outname,'recreate')
                title = 'pic_0'
                h2 = ROOT.TH2F(title,title,width,0,width,height,0,height)
                h2.GetXaxis().SetTitle('x')
                h2.GetYaxis().SetTitle('y')
                for iy in range(height):
                    for ix in range(width):
                        h2.SetBinContent(ix+1,iy+1,data[iy,ix])
                h2.Write()
                rf.Close()

                # Plot the 2D histogram
                plt.imshow(data, cmap='viridis', interpolation='nearest')
                plt.colorbar(label='Value')
                plt.title('2D Histogram')
                plt.xlabel('X-axis')
                plt.ylabel('Y-axis')
                plt.show()

            else:
                print("No numbers found in the file.")

    except FileNotFoundError:
        print(f"The file '{file_path}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <file_path>")
    else:
        file_path = sys.argv[1]
        process_file(file_path)
