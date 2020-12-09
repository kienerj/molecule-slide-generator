import unittest
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from molecule_slide_generator import *
from PIL import Image
import imagehash


class ImageCreationTest(unittest.TestCase):
    """
    Use ImageHash library to compare newly generated image to a pre-generated image
    This gives no guarantee the image is 100% correct but prevents from manual testing for every change
    """

    def test_default_image(self):

        sg = SlideGenerator(number_of_properties=2)
        png = sg.generate_slide(self.mols, self.all_props, 'test_default.png')
        hash = imagehash.dhash(Image.open('test_default.png'), hash_size=16)
        self.assertEqual(self.default_hash, hash, msg="Images hashes do not match.")

        # Check properties
        num_properties = 2
        # pickle, molfile and smiles
        num_mol_formats = 3
        num_mols = 21
        # +1 for num_properties entry in metadata
        num_keys = (num_properties + num_mol_formats) * num_mols + 1
        props = Chem.MetadataFromPNGString(png)
        self.assertEqual(num_keys, len(props.keys()), msg="Number of key in metadata dict does not match. Expected {}, "
                                                          "got {}.".format(num_keys, len(props.keys())))
        self.check_extraction(png, num_mols)

    def test_example_image(self):
        sg = SlideGenerator(mols_per_row=3, rows=3, font_size=18, font='comicbd', number_of_properties=2,
                            slide_width=800, slide_height=600)
        png = sg.generate_slide(self.mols, self.all_props, 'test_example.png')
        hash = imagehash.dhash(Image.open('test_example.png'), hash_size=16)
        self.assertEqual(self.example_hash, hash, msg="Images hashes do not match.")

        # Check properties
        num_properties = 2
        # pickle, molfile and smiles
        num_mol_formats = 3
        num_mols = 9
        # +1 for num_properties entry in metadata
        num_keys = (num_properties + num_mol_formats) * num_mols + 1
        props = Chem.MetadataFromPNGString(png)
        self.assertEqual(num_keys, len(props.keys()), msg="Number of key in metadat dict does not match. Expected {}, "
                                                          "got {}.".format(num_keys, len(props.keys())))

        self.check_extraction(png, num_mols)

    def check_extraction(self, png, num_mols):

        e = DataExtractor()
        df = e.extract_single(png)
        self.assertEqual(num_mols, len(df),
                         msg="Extracted dataframe has wrong length. Expected {}, got {}.".format(num_mols, len(df)))

        name0 = df.iloc[0, 0]
        self.assertEqual('Adinazolam', name0,
                         msg="First molecule has wrong name. Expected Adinazolam, got {}.".format(name0))
        name5 = df.iloc[5, 0]
        self.assertEqual('Clonazepam', name5,
                         msg="First molecule has wrong name. Expected Clonazepam, got {}.".format(name5))

        mol0 = df.iloc[0,2]
        smi = Chem.MolToSmiles(mol0)
        self.assertEqual('CN(C)Cc1nnc2CN=C(c3ccccc3)c3c(ccc(Cl)c3)-n21', smi, msg="First molecule does not match.")

    def setUp(self):
        suppl = Chem.SDMolSupplier('bzr21.sdf')
        self.mols = [x for x in suppl]
        self.all_props = []
        for mol in self.mols:
            AllChem.Compute2DCoords(mol)
            a = float(mol.GetProp('ACTIVITY'))
            if a > 8.0:
                color = '#3A662F'
            elif a < 8.0 and a > 7:
                color = '#e8860e'
            else:
                color = '#b52009'
            props = [TextProperty('Name', mol.GetProp('_Name')), TextProperty('Activity', a, color=color)]
            self.all_props.append(props)

        self.default_hash = imagehash.dhash(Image.open('../images/test_default.png'), hash_size=16)
        self.example_hash = imagehash.dhash(Image.open('../images/example_slide_bzr.png'), hash_size=16)


if __name__ == '__main__':
    unittest.main()
