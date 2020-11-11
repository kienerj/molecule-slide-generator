from rdkit import Chem
from rdkit.Chem import PandasTools
from io import BytesIO
import io
from PIL import Image
import re
import pandas as pd
import io


class DataExtractor(object):
    """
    Extracts the molecules and properties out of the image metadata into a pandas dataframe

    The extraction can be done on a single image (extract_single) or repeatedly called in a loop (extract). Properties
    not available for a specific molecule will be set to `None`.

    If a molecule doesn't contain a molecule, it is silently ignored.
    """

    def __init__(self):

        self.reference = []
        self.molecules = []
        self.columns = []
        self.data = {}
        self.image_number = 0

    def extract_single(self, image):
        """
        Extracts molecules and properties from a single png image

        :param image: image bytes, file path or binary stream
        :return: DataFrame with the molecules and properties
        """

        self.extract(image)
        return self.get_data()

    def extract(self, image):
        """
        Extracts molecules and properties from a png image.

        This method can be called repeatedly in a loop and when all images are extracted,
        the data can be fetched by a call to `get_data()`.

        :param image: image bytes, file path or binary stream
        """

        if isinstance(image, bytes):
            img = image
        elif isinstance(image, str):
            with open(image, 'rb') as f:
                img = f.read()
        elif isinstance(image, io.BufferedIOBase):
            img = image.read()
        else:
            raise ValueError("Can't determine data type of image. Supported are bytes, str or a binary stream.")

        stream = BytesIO(img)
        pil_image = Image.open(stream)
        pil_image.load()
        is_molimage = False
        for key, value in pil_image.info.items():
            if key.startswith('rdkit'):
                is_molimage = True
                break
        if not is_molimage:
            # no molecule in image, do not process
            return

        mols = Chem.MolsFromPNGString(img)
        nprops = Chem.MetadataFromPNGString(img)

        if 'numProperties' in nprops:
            # remove structure data properties which all contain 'rdkit'
            rx = re.compile(r'rdkit')
            props = {key: nprops[key] for key in nprops if not rx.search(key)}
            num_props = int(props.pop('numProperties').decode('utf-8'))

            prop_names = list(props)
            for idx in range(num_props):
                name = prop_names[idx]
                if name not in self.columns:
                    self.columns.append(name)
                    self.data[name] = []  # initialize data
                    if len(self.molecules) > 0:
                        # not first image, fill list with none
                        self.data[name] = [None] * len(self.molecules)
        else:
            props = {}
        for idx in range(len(mols)):
            for prop in self.columns:
                if idx > 0:
                    key = prop + str(idx)
                else:
                    key = prop
                if key in props:
                    self.data[prop].append(props[key].decode('utf-8'))
                else:
                    self.data[prop].append(None)

        self.molecules.extend(mols)  # extend here to simplify above checks
        self.reference.extend([self.image_number] * len(mols))
        self.image_number += 1

    def get_data(self):

        self.data['Molecule'] = self.molecules
        self.data['Reference'] = self.reference
        df = pd.DataFrame(data=self.data)
        PandasTools.ChangeMoleculeRendering(df)
        self.__init__() #reinitalize
        return df

