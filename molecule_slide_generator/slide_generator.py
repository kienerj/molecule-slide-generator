from PIL import Image, ImageDraw, ImageFont
from PIL.PngImagePlugin import PngInfo
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from io import BytesIO
import math
import numpy as np
import rdkit
import platform
import os
import pathlib


class SlideGenerator(object):
    """Molecule Overview Slide Generator

    This class creates png images meant to be used in presentations that need simple overviews of molecules.
    Each slide (image) contains a certain amount of molecules plus additional properties of interest (text) below the
    molecule like molecule name, internal registration number and activity data.

    The generator also stores the molecule info and the properties inside the pngs metadata and they can therefore be
    extracted again at a later time if so desired. (The reason png is used is because MS Office products strip out
    metadata out of the svg but not out of pngs.)

    """

    # Font folder on Ubuntu-type distros
    LINUX_TRUETYPE_FONTS_FOLDER = '/usr/share/fonts/truetype'

    def __init__(self, mols_per_row=7, rows=3, font_size=16, font='Calibrib', number_of_properties=4, bond_length=20,
                 bond_width=2, slide_width=1440, slide_height=607, dpi = (120,120)):
        """
        The parameters defined are used to calculate the size that each molecule and the properties part can take up.

        The defaults result in a well looking overview of molecules for putting into a PowerPoint presentation that
        has a 16:9 aspect ratio and uses a Title->Subtitle slide template. Default slide_width and slide_height were
        determined manually and optimal values will depend on your exact templates. PowerPoint respects the DPI metadata
        of the generated png image and scales it accordingly. So for different templates/use cases you will need to do
        some manual experiments to figure out optimal settings.

        'font' must either be a font name of a installed font or the full path to a true type font file (.ttf).
        For Windows the named fonts must be in 'C:/Windows/Fonts/' directory and for Linux in '/usr/share/fonts/truetype'.
        (Linux font folder is for Ubuntu. Don't know if it is the same for other distros).
        The font is used for the atom labels and the properties.

        'number_of_properties' defines how many additional data you want to show below the molecules. Each property takes
        up 1 line of text. Example: pIC50: 6.2. This is meant for very small amounts of text/data. If text is too long,
        it is simply truncated automatically due to image width limit. It's up to the caller to adjust the data if hard
        truncation isn't acceptable. The approximate number of characters per line can be deduced from the instance
        attribute 'num_chars_per_line'. Exact amount depends on the exact text.

        :param mols_per_row: number of molecules per row. default:7
        :param rows: number of rows. default: 3
        :param font_size: size of the properties font. default:16
        :param font: name of the font to use or path to a .ttf font file
        :param number_of_properties: number of properties to display below the molecule
        """

        root, ext = os.path.splitext(font)
        platform_system = platform.system()
        if ext != '' and ext != '.ttf':
            raise ValueError("'font' must be a path to an existing ttf file. However {} was provided as font.".format(font))
        elif ext == '.ttf':
            self.font_path = font
        elif platform_system == 'Windows':
            self.font_path = 'C:/Windows/Fonts/' + font + '.ttf'
        elif platform_system == 'Linux':
            # For Windows fonts, requires MS Core fonts installed
            # sudo apt install ttf-mscorefonts-installer
            font_file_name = font + '.ttf'
            for root, dirs, files in os.walk(SlideGenerator.LINUX_TRUETYPE_FONTS_FOLDER):
                if font_file_name in files:
                    self.font_path = os.path.join(root, font_file_name)
                    break
        else:
            raise ValueError ("'font' was set to name {} but the location of the font folder is unknown for platform "
                              "{}. Please use full path to font file.".format(font, platform_system))

        p = pathlib.Path(self.font_path)
        if not p.is_file():
            raise ValueError("Desired font file '{}' does not exist.".format({self.font_path}))

        self.font_size = font_size
        self.bond_width = bond_width
        self.mols_per_row = mols_per_row
        self.mols_per_column = rows
        self.max_mols = mols_per_row * rows
        self.number_of_properties = number_of_properties
        self.slide_width = slide_width
        self.slide_height = slide_height
        self.dpi = dpi
        self.image_width = slide_width // mols_per_row
        self.row_height = slide_height // rows
        self.font = ImageFont.truetype(self.font_path, size=font_size)
        self.line_height = self.font.getsize('hg')[1] + 1
        self.text_image_height = math.ceil(self.line_height * number_of_properties)
        self.molecule_image_height = self.row_height - self.text_image_height
        # estimate max number of characters per line
        avg_char_width = self.font.getsize('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ123456789-, ')[0] / 64
        self.num_chars_per_line = self.image_width / avg_char_width
        self.bond_length = bond_length

    def generate_slide(self, mols, properties, out_path=None):
        """
        Generates an image with all molecules and their properties below them.

        :param mols: list of rdkit molecules
        :param properties: list of list of TextProperty
        :param out_path: optional file to save generated images
        :return: bytes object containing the png image data say for display in a notebook
        """

        mols, slide, png_info = self._prepare_slide(mols, properties)

        for index, mol in enumerate(mols):

            mol_image = self._draw_mol(mol)

            self._add_to_slide(slide, properties, mol_image, index, png_info)

        return self._finalize_slide(slide, png_info, out_path)

    def generate_slide_from_images(self, mols, properties, images, out_path=None):
        """
        Generates an image with all molecules and their properties below them.

        Images is a list containing the molecules images as PIL images to use in the slide. They must be in same order as mols.
        This allows each molecule image to be generated differently like substructure highlighting, annotations and so
        forth. So many options they can't be all part of a methods argument.

        :param mols: list of rdkit molecules
        :param properties: list of list of TextProperty
        :param images: list of images to use for each molecule
        :param out_path: optional file to save generated images
        :return: bytes object containing the png image data say for display in a notebook
        """

        if images is None or len(images) != len(mols):
            raise ValueError("List of images is of different size than list of molecules.")

        mols, slide, png_info = self._prepare_slide(mols, properties)

        for index, mol in enumerate(mols):

            mol_image = images[index]
            if mol_image.width != self.image_width:
                raise ValueError("Image at index {} has wrong width. {} but expected {}."
                                 .format(index, mol_image.width, self.image_width))
            if mol_image.height != self.molecule_image_height:
                raise ValueError("Image at index {} has wrong height. {} but expected {}."
                                 .format(index, mol_image.height, self.molecule_image_height))

            self._add_to_slide(slide, properties, mol_image, index, png_info)

        return self._finalize_slide(slide, png_info, out_path)

    def _prepare_slide(self, mols, properties):

        if mols is None:
            raise ValueError('Expected a list of rdkit molecule instances but got \'None\'')

        # cut-off mols + properties silently. I think this is better than raising a ValueError
        mols = mols[:self.max_mols]
        properties = properties[:self.max_mols]
        if self.number_of_properties > 0 and len(mols) != len(properties):
            raise ValueError('Number of molecules must match number of properties.')

        if len(mols) == self.max_mols:
            slide = Image.new('RGBA', [self.slide_width, self.slide_height], (255, 255, 255, 0))
        # less mols -> smaller image
        elif len(mols) < self.mols_per_row:
            slide = Image.new('RGBA', [len(mols) * self.image_width, self.row_height], (255, 255, 255, 0))
        else:
            num_rows = (len(mols) // self.mols_per_row) + 1
            slide = Image.new('RGBA', [self.slide_width, self.row_height * num_rows], (255, 255, 255, 0))

        png_info = PngInfo()
        png_info.add_text('numProperties', str(self.number_of_properties))

        return mols, slide, png_info

    def _add_to_slide(self, slide, properties, mol_image, index, png_info):

        if self.number_of_properties > 0:
            text_image = self._draw_text(properties[index])

            combined = Image.new('RGBA', [self.image_width, self.row_height], (255, 255, 255, 0))
            combined.paste(mol_image)
            combined.paste(text_image, (0, self.molecule_image_height))
        else:
            combined = mol_image

        row = index // self.mols_per_row
        column = index % self.mols_per_row
        position = (column * self.image_width, row * self.row_height)

        slide.paste(combined, position)

        if index > 0:
            key_index = str(index)
        else:
            key_index = ''

        png_info.add_text('rdkitPKL{} rdkit {}'.format(key_index, rdkit.__version__),
                          mol_image.info['rdkitPKL rdkit {}'.format(rdkit.__version__)])
        png_info.add_text('MOL{} rdkit {}'.format(key_index, rdkit.__version__),
                          mol_image.info['MOL rdkit {}'.format(rdkit.__version__)])
        png_info.add_text('SMILES{} rdkit {}'.format(key_index, rdkit.__version__),
                          mol_image.info['SMILES rdkit {}'.format(rdkit.__version__)])

        if self.number_of_properties > 0:
            mol_properties = properties[index]
            mol_properties = mol_properties[:self.number_of_properties]
            for prop in mol_properties:
                png_info.add_text('{}{}'.format(prop.name, key_index), str(prop.value))
                # print("{}: {}".format(prop.name, str(prop.value)))

    def _finalize_slide(self, slide, png_info, out_path):

        if out_path is not None:
            slide.save(out_path, format='PNG', dpi=self.dpi, pnginfo=png_info)

        img_byte_arr = BytesIO()
        slide.save(img_byte_arr, format='PNG', dpi=self.dpi, pnginfo=png_info)
        img_byte_arr = img_byte_arr.getvalue()
        return img_byte_arr

    def _draw_text(self, mol_properties):

        mol_properties = mol_properties[:self.number_of_properties]
        txt = Image.new('RGBA', (self.image_width, self.text_image_height), (255, 255, 255, 0))
        d = ImageDraw.Draw(txt)
        position = (0, 0)

        for prop in mol_properties:
            d.text(position, prop.get_display_value(), font=self.font, fill=prop.color)
            position = (0, position[1] + self.line_height)

        return txt

    def _draw_mol(self, mol, kekulize=True, add_chiral_hs=True):

        if mol is None:
            raise ValueError('Expected a rdkit molecule instance but got \'None\'')

        mol = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=kekulize, addChiralHs=add_chiral_hs)
        mol = SlideGenerator._scale_bond_length(mol)
        drawer = rdMolDraw2D.MolDraw2DCairo(self.image_width, self.molecule_image_height)
        drawer.SetFontSize(self.font_size)
        opts = drawer.drawOptions()
        opts.clearBackground = False
        opts.bondLineWidth = self.bond_width
        opts.fixedBondLength = self.bond_length
        opts.minFontSize = self.font_size
        opts.maxFontSize = self.font_size
        opts.fontFile = self.font_path

        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        png = drawer.GetDrawingText()
        mol_image = Image.open(BytesIO(png))
        mol_image.load()
        return mol_image

    @staticmethod
    def _scale_bond_length(mol):

        default_bond_length = 1.5  # rdkit default bond length

        bonds = mol.GetBonds()

        if len(bonds) == 0:
            return mol

        total = 0
        for bond in bonds:

            ai = bond.GetBeginAtomIdx()
            aj = bond.GetEndAtomIdx()
            bl = AllChem.GetBondLength(mol.GetConformer(), ai, aj)
            if not math.isnan(bl):
                total += bl

        avg_bl = (total / len(bonds))

        if avg_bl > 0.0:

            bl_ratio = default_bond_length / avg_bl

            tm = np.zeros((4, 4), np.double)
            for i in range(3):
                tm[i, i] = bl_ratio

            AllChem.TransformMol(mol, tm)
            return mol
        else:
            # in some cases avg_bl is zero -> recreate 2D coords and try again
            AllChem.Compute2DCoords(mol)
            return SlideGenerator._scale_bond_length(mol)


class TextProperty(object):

    def __init__(self, name, value, show_name=False, color=(0, 0, 0, 255)):
        """
        Sets options of the data below the molecule is displayed.

        Color is either a HEX value or a RGBA 4-tuple.

        :param name: display name of the property
        :param value: the value of the property
        :param show_name: if the name of the property should be displayed or not
        :param color: color of the text for this property. Default is black.
        """
        self.name = name
        self.value = value
        self.show_name = show_name

        if isinstance(color, tuple):
            if len(color) == 4:
                self.color = color
            else:
                raise ValueError('Expected a RGBA color tuple of 4 values. Got tuple with {} values'.format(len(color)))
        elif isinstance(color, str) and color[0] == '#':
            self.color = TextProperty.hex_to_rgb(color)
        else:
            raise ValueError('Expected a hex color string or RGBA 4-tuple but got {}.'.format(color))

    def get_display_value(self):

        if self.show_name:
            return self.name + ': ' + str(self.value)
        else:
            return str(self.value)

    @staticmethod
    def hex_to_rgb(hex_code):
        # copied from stackoverflow
        hex_code = hex_code.lstrip('#').upper()
        rgb = tuple(int(hex_code[i:i + 2], 16) for i in (0, 2, 4))
        return rgb + (255,)  # RGBA
