# RDKit Molecule Slide Generator

RDKit Molecule Slide Generator is a tool for generating an overview image of molecules and their properties for usage in reports or presentations. It's easier to show than explain.

![example_slide_bzr](images/example_slide_bzr.png)

This examples uses the first 9 molecules of the `bzr.sdf` in the RDKit distributions data directory and showcases a number of features.

### Features

- Any true type font (.ttf) can be used for the atom labels and text. Comic Sans MS in this case.
- Set font size of atom labels and text
- The text can be colored separately for each molecule
- Set bond length of molecule (`fixedBondLength`drawing option)
- The property can be displayed with the properties name instead of just the value as above
- Define size of slide (width, height) and number of rows and columns 
- Molecule data and properties are stored in images metadata and can be extracted again

### Requirements

- RDKit 2020.09.1 or higher (due to new PNG Metadata feature used here)
- numpy
- pillow

### Example Code

This is the code used to generate above image on Windows. For Linux you would have to use a different font or install MS Core fonts.

```python
suppl = Chem.SDMolSupplier('data/bzr.sdf')
# first 9 mols
mols = [x for x in suppl][:9]

# Create properties with coloring
all_props = []
for mol in mols:    
    AllChem.Compute2DCoords(mol) # remove 3D coords
    a = float(mol.GetProp('ACTIVITY'))
    if a > 8.0:
        color = '#3A662F'
    elif a < 8.0 and a > 7:
        color = '#e8860e'
    else:
        color = '#b52009'
    props = [TextProperty('Name', mol.GetProp('_Name')), 
             TextProperty('Activity', a, color=color)]
    all_props.append(props)
    
# Generate the slide
sg = SlideGenerator(mols_per_row=3, rows=3, font_size=18, font='comicbd', 
                    number_of_properties=2, slide_width=800, slide_height=600)
png = sg.generate_slide(mols, all_props, 'example_slide.png')
```

For font you can either use a named font or a full path to a `.ttf` file. For named fonts the font must be in `C:/Windows/Fonts/` directory for Windows and for Linux in `/usr/share/fonts/truetype`.

To read out the molecules, use the according RDKit functionality:

```python
mols = Chem.MolsFromPNGString(png)
Draw.MolsToGridImage(mols,molsPerRow=3)
```

where `png`is a bytes object containing the png data.

For properties you can see the available data also by using RDKit:

```python
nprops = Chem.MetadataFromPNGString(png)
nprops.keys()
dict_keys(['rdkitPKL rdkit 2020.09.1', 'MOL rdkit 2020.09.1', 'SMILES rdkit 2020.09.1', 'Name', 
           'Activity', 'rdkitPKL1 rdkit 2020.09.1', 'MOL1 rdkit 2020.09.1', 'SMILES1 rdkit 2020.09.1', 
           'Name1', 'Activity1', 'rdkitPKL2 rdkit 2020.09.1', 'MOL2 rdkit 2020.09.1', 
           'SMILES2 rdkit 2020.09.1', 'Name2', 'Activity2', 'rdkitPKL3 rdkit 2020.09.1', 
           'MOL3 rdkit 2020.09.1', 'SMILES3 rdkit 2020.09.1', 'Name3', 'Activity3', 
           'rdkitPKL4 rdkit 2020.09.1', 'MOL4 rdkit 2020.09.1', 'SMILES4 rdkit 2020.09.1', 'Name4', 
           'Activity4', 'rdkitPKL5 rdkit 2020.09.1', 'MOL5 rdkit 2020.09.1', 'SMILES5 rdkit 2020.09.1', 
           'Name5', 'Activity5', 'rdkitPKL6 rdkit 2020.09.1', 'MOL6 rdkit 2020.09.1', 
           'SMILES6 rdkit 2020.09.1', 'Name6', 'Activity6', 'rdkitPKL7 rdkit 2020.09.1', 
           'MOL7 rdkit 2020.09.1', 'SMILES7 rdkit 2020.09.1', 'Name7', 'Activity7', 
           'rdkitPKL8 rdkit 2020.09.1', 'MOL8 rdkit 2020.09.1', 'SMILES8 rdkit 2020.09.1', 'Name8', 
           'Activity8'])
```

The system used is that first molecules data is in key with the properties name, eg. "Activity" in the example. For each subsequent molecule the index of the molecule in the file is appended. So to get the activity value for the second molecule in the image you would use the key "Activity1". Data is returned as bytes and must be converted accordingly.

```python
float(nprops['Activity1'].decode('utf-8'))
7.7
```

### Miscellaneous

The default parameters of the slide generator work well with small molecules and a PowerPoint template with an aspect ratio of 16:9 that has a title and subtitle text box, the specific setup I have to use. You will need to play with the parameters to find an good solution for your specific needs.

### Exporting the data from MS Word or PowerPoint

If you insert the image into an MS Office document, it can get problematic to extract the original image with all the metadata. However looking at the document contents the original images are present and can be extracted again. For extracting all the molecules in an MS Office document you can use the KNIME component [Extract RDKit Molecules From Office](https://hub.knime.com/kienerj/spaces/Public/latest/Extract%20RDKit%20Molecules%20From%20Office).