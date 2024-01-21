import marimo

__generated_with = "0.1.79"
app = marimo.App()


@app.cell
def __():
    # Read the final geometries and calculate the volume of the cells
    from ase.io import read
    from ase.visualize import view
    return read, view


@app.cell
def __():
    # Make a list of filenames in the current directory
    import os
    filenames = [f for f in os.listdir('.') if f.endswith('.pdb')]
    return filenames, os


@app.cell
def __(filenames):
    # Sort filenames in natural order
    import re
    filenames.sort(key=lambda f: int(re.sub('\D', '', f)))
    return re,


@app.cell
def __(filenames):
    filenames
    return


@app.cell
def __(filenames, read, view):
    view([read(f) for f in filenames])
    return


@app.cell
def __(filenames):
    # Extract the pressure from each filename
    pressures = [(f.split('_')[1]) for f in filenames]
    # Strip the file extension from pressures
    pressures = [float(p[1:-7]) for p in pressures]
    pressures
    return pressures,


@app.cell
def __(filenames, read):
    # Define the array of volumes
    volumes = []
    # Read geometries for each file
    for i, f in enumerate(filenames):
        # Read the geometry
        atoms = read(f)
        # Calculate the volume
        volume = atoms.get_volume()
        # Append the volume to the list
        volumes.append(volume)
    return atoms, f, i, volume, volumes


@app.cell
def __(pressures, volumes):
    # Plot the volume as a function of pressure
    import matplotlib.pyplot as plt
    plt.plot(pressures, volumes, 'o')
    plt.xlabel('Pressure (GPa)')
    plt.ylabel('Volume ($\AA^3$)')
    plt.show()
    return plt,


@app.cell
def __(plt):
    # Fetch experimental data on the volume of ice as a function of pressure
    import urllib.request
    url = 'http://www1.lsbu.ac.uk/water/water_phase_diagram.html'
    with urllib.request.urlopen(url) as response:
        html = response.read()
        # Parse the HTML
        from bs4 import BeautifulSoup
        soup = BeautifulSoup(html, 'html.parser')
        # Extract the table
        table = soup.find_all('table')[0]
        # Extract the rows
        rows = table.find_all('tr')
        # Extract the data from each row
        data = []
        for row in rows:
            cols = row.find_all('td')
            cols = [c.text.strip() for c in cols]
            data.append(cols)
            # Convert the data to a Pandas DataFrame
            import pandas as pd
            df = pd.DataFrame(data)
            # Set the column names
            df.columns = df.iloc[0]
            # Drop the first row
            df = df.reindex(df.index.drop(0))
            # Convert the data to floats
            df = df.astype(float)
            # Plot the experimental data
            plt.plot(df['Pressure (GPa)'], df['Volume ($\AA^3$)'], 'o')
            plt.xlabel('Pressure (GPa)')
            plt.ylabel('Volume ($\AA^3$)')
            plt.show()
            
    return (
        BeautifulSoup,
        cols,
        data,
        df,
        html,
        pd,
        response,
        row,
        rows,
        soup,
        table,
        url,
        urllib,
    )


if __name__ == "__main__":
    app.run()
