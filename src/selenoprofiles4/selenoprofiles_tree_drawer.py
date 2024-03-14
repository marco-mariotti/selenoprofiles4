#!/home/mmariotti/miniconda3/envs/selenoprofiles/bin/python -u
from string import *
import sys
from subprocess import *
from .MMlib3 import *
from .selenoprofiles4 import (
    p2ghit,
    blasthit,
    exoneratehit,
    genewisehit,
    notracebackException,
)

try:
    from ete3 import *
    from PyQt5 import QtCore, QtGui, QtWidgets

except:
    # printerr
    raise notracebackException(
        "ERROR ete3 and PyQt5 must be installed to use selenoprofiles drawer! Try this command:\nconda install -c etetoolkit ete3\n\nor follow instructions at http://etetoolkit.org/download/"
    )

try:
    from PyQt5.QtWidgets import (
        QGraphicsRectItem,
        QGraphicsLineItem,
        QGraphicsSimpleTextItem,
    )

except:
    from PyQt5.QtGui import (
        QGraphicsRectItem,
        QGraphicsLineItem,
        QGraphicsSimpleTextItem,
    )

from PyQt5.QtGui import QBrush, QPen, QColor

help_msg = """selenoprofiles drawer: utility to draw the distribution of selenoprofiles results on a species tree. 

Usage:  selenoprofiles drawer -t tree.nw -i family1.ali  [family2.ali .. familyM.ali] [options]

Input alignments must be .ali files as produced by selenoprofiles join.
In alignments, only selenoprofiles results are considered. Species names and intron coordinates are derived the title.
Every species found in ali files must be present in the tree, or the program crashes (unless option -g is active).
The tree file can be omitted to draw an unstructured tree.

The tree can contain also species without result. These nodes are not drawn, unless option -e is active.

In the output, results are colored by their selenoprofiles label, e.g. "selenocysteine" results are green. 
Use the -colors option to display the correspondance between labels and colors.

The program opens an interactive ETE3 graphics unles option -out is active.
If you're on a remote server, use X11 tunneling in your ssh connection.

* I/O options (+ means an argument is required):
-i +      space-separated ali file(s) from selenoprofiles join
-t +      species tree in newick format
-m        the input tree has masked species names; e.g. Xenopus (Silurana) tropicalis is coded as Xenopus_{ch40}Silurana{ch41}_tropicalis
-g        do not crash if some ali results are missing from the tree. These are ignored, and a warning is printed for each
-out +    save image in the output file provided as argument (allowed extensions: pdf, png, svg). The interactive environment is not opened

* Tree graphics options:
-e        empty nodes (without any assigned result) are also drawed
-o  +     outgroup species. If not provided, the midpoint species is used. Use -o None to not set it
-img_w +  image width, used for -out option
-img_h +  image height, used for -out option
-C        draw as a circular tree
-F        print family names in each species line, not just as header of their columns 
-a        suppress normal output; colored numbers are used to summarize the number of hits per label for each family
"""


help_msg_full = (
    help_msg
    + """  # if option -a is active:
  -c                compress the number of boxes by summing up all labels that are not "explicit labels"
  -explicit_labels  define explicit labes for option -c; format: -explicit_labels label1,label2,label3

* Graphics for each gene element, applicable if option -a is NOT active:
-W  +     gene brick width in pixels (default: 150)
-H  +     gene brick height in pixels (default: 40)
-T        don't display text on gene bricks (chromosome names and positions)
-I        don't display intron markers 
-no_id    don't display numeric id on the left
-f +      display colored vertical lines for features indicated in fasta headers. Argument form: feat1,feat2,feat3
          Each feature is searched in fasta headers, looking for something like (in case of feat1):  " feat1:1,33,41 "
          These numbers are understood as positions relative to the protein sequence. 
          Optionally you can define line colors adding a RGB code after the feature name, e.g  feat1#33FF44,feat2#EE0000
-add  +   provide an tab delimited file to display additional features for the results. This file must have an entry per line, like protein_idTABtext 
-sp_add + add annotation for each species from a tab delimited file, format: 'species TAB annotation'. 
          The annotations are appended next to the species name (aligned). 
          The special annotation COLOR#xxxxx can be provided to change the color of the species name to xxxxxx (RGB, hex)

* Text size:   
-tsize +  specify size of text for titles
-ssize +  specify size of text for species 
-bsize +  specify size of body text
-nsize +  specify size of numbers in the colored boxes (when option -a is active)
NOTE: if a text doesn't not fit its dedicated space, its point size is reduced until it does.

* Other options:
-prompt       open interactive prompt at the end
-print_opt    print currently active options
-h OR --help  print this help and exit
"""
)

help_msg += "\nThere are other options to set graphical parameters. Display them with: selenoprofiles drawer -h full"

""" -common   convert species names from scientific to common (only for species for which it's hardcoded -- see common_names hash in the script code)
"""


def set_selenoprofiles_tree_drawer_var(varname, value):
    globals()[varname] = value


command_line_synonyms = {}

def_opt = {  #'temp':'/home/mmariotti/temp', 'common':0,
    "colors":False,
    "no_id": False,
    "add": 0,
    "sp_add": 0,
    "t": "",  # "/users/rg/mmariotti/Genomes/.tree",
    "prompt": False,
    "m": False,
    #'s':'',
    "out": "",
    "a": False,
    "C": False,
    "c": False,
    "explicit_labels": "selenocysteine,cysteine,homologue",
    "I": False,
    "F": False,
    "e": False,
    "g": False,
    "T": False,
    "f": False,
    "o": "",
    #'O':'',
    "tsize": 20,
    "ssize": 20,
    "bsize": 8,
    "nsize": 14,
    "margin_boxes": 10,
    "W": 150,
    "H": 40,
    #'v':0, 'Q':0, 'ali':0, 'abs_pos':0,
    "intron_color": "#FFFFFF",
    "img_h": -1,
    "img_w": -1,
    "cmd": "drawer",
    "i": [],
    "filter":False
}


#########################################################
###### start main program function

global label_to_color
color_tuples=[
    ("selenocysteine",  "#89C900",  "green"),
    ("cysteine",  "#F40000", "red" ),
    ("arginine",  "#F15BA6", "pink" ),
    ("threonine",  "#8B1B8D", "dark purple" ),
    ("pseudo",  "#545454", "gray" ),
    ("uga_containing",  "#EEA347", "orange" ),
    ("unaligned",  "#979796", "grey" ),
    ("gapped",  "#979796", "grey" ),
    ("homologue",  "#E1D100", "yellow" ),
    ("unknown",  "#979796", "light grey" ),
    ("serine",  "#D09E5F", "light brown" ),
    ("readthrough",  "#0E6DBB", "blue" ),
    ("tRNA",  "#C5D51B", "yellow" ),
    ("glycine",  "#FFB20B", "bright orange" ),
    ("leucine",  "#834D1A", "brown" ),
    ("tRNA?",  "#7A850D", "dark yellow" ),
    ("Missing","#FFFF66","light salmon")
    ]

label_to_color={x[0]:x[1] for x in color_tuples}


global ordered_labels
ordered_labels = [
    "selenocysteine",
    "cysteine",
    "arginine",
    "threonine",
    "uga_containing",
    "unaligned",
    "gapped",
    "alanine",
    "asparagine",
    "aspartic_acid",
    "glutamic_acid",
    "glutamine",
    "glycine",
    "histidine",
    "isoleucine",
    "leucine",
    "lysine",
    "methionine",
    "phenylalanine",
    "proline",
    "serine",
    "tryptophan",
    "tyrosine",
    "valine",
    "homologue",
    "pseudo",
    "readthrough",
    "tRNA",
    "tRNA?",
    "Missing"
]

common_names = {
    "Ornithorhynchus_anatinus": "Platypus",
    "scientific_name": "common name",
    "Oryctolagus_cuniculus": "Rabbit",
    "Dipodomys_ordii": "Kangaroo_rat",
    "Sus_scrofa": "Pig",
    "Pongo_pygmaeus": "Orangutan",
    "Sorex_araneus": "Shrew",
    "Myotis_lucifugus": "Microbat",
    "Oryzias_latipes": "Medaka",
    "Monodelphis_domestica": "Opossum",
    "Anolis_carolinensis": "Lizard",
    "Callorhinchus_milii": "Elephant_shark",
    "Microcebus_murinus": "Mouse_lemur",
    "Lama_pacos": "Llama",
    "Taeniopygia_guttata": "Finch",
    "Pan_troglodytes": "Chimpanzee",
    "Gasterosteus_aculeatus": "Stickleback",
    "Homo_sapiens": "Human",
    "Tupaia_belangeri": "Tree_shrew",
    "Dasypus_novemcinctus": "Armadillo",
    "Macaca_mulatta": "Macaque",
    "Otolemur_garnettii": "Galago",
    "Spermophilus_tridecemlineatus": "Squirrel",
    "Rattus_norvegicus": "Rat",
    "Macropus_eugenii": "Wallaby",
    "Xenopus_{ch40}Silurana{ch41}_tropicalis": "Frog",
    "Gallus_gallus": "Chicken",
    "Erinaceus_europaeus": "Hedgehog",
    "Gorilla_gorilla": "Gorilla",
    "Mus_musculus": "Mouse",
    "Choloepus_hoffmanni": "Sloth",
    "Tursiops_truncatus": "Dolphin",
    "Takifugu_rubripes": "Fugu",
    "Felis_catus": "Cat",
    "Callithrix_jacchus": "Marmoset",
    "Bos_taurus": "Cow",
    "Equus_caballus": "Horse",
    "Canis_lupus_familiaris": "Dog",
    "Pteropus_vampyrus": "Macrobat",
    "Danio_rerio": "Zebrafish",
    "Procavia_capensis": "Hyrax",
    "Loxodonta_africana": "Elephant",
    "Cavia_porcellus": "Guinea_pig",
    "Tetraodon_nigroviridis": "Pufferfish",
    "Tarsius_syrichta": "Tarsier",
    "Branchiostoma_floridae": "Lancelet",
    "Drosophila_melanogaster": "Fruit_fly",
    "Anopheles_gambiae": "Malaria_mosquito",
    "Aedes_aegypti": "Yellow_fever_mosquito",
    "Acyrthosiphon_pisum": "Pea_aphid",
    "Pediculus_humanus": "Louse",
    "Tribolium_castaneum": "Beetle",
    "Apis_mellifera": "Honey_bee",
    "Nasonia_vitripennis": "Jewel_wasp",
}

global warnings_reduced_text
warnings_reduced_text = {}  # keys: categories


def reduce_font_if_necessary(simpleTextItem, width=-1, height=-1, category="unknown"):
    psize = 99  # not used, this is just for the first time last condition of while is checked
    reducing_of_points = 0
    while (
        (width > 0 and simpleTextItem.boundingRect().width() > width)
        or (height > 0 and simpleTextItem.boundingRect().height() > height)
        and psize > 0
    ):
        font = simpleTextItem.font()
        psize = font.pointSize()
        font.setPointSize(psize - 1)
        simpleTextItem.setFont(font)
        reducing_of_points += 1
    if reducing_of_points:
        if category not in warnings_reduced_text:
            warnings_reduced_text[category] = 0
        warnings_reduced_text[category] += 1
        # printerr('Reduced point size for text "'+simpleTextItem.text()+'" in category: '+category, 1)


class NumberedBoxFace(faces.Face):
    """This face displays numbers in colored boxes. Initialise it with a list like [ [number, '#hexcodeforcolor'], [], [number, '#hexcodeforcolor']  ...  ]
  Each element of this list is a two-element list with the number to write and the corresponding color in hex code. Empty elements are display like white boxes (nothing is drawn in this position) this is useful to keep boxes in different lines still aligned 
 """

    def __init__(self, input_list):
        faces.Face.__init__(self)
        self.type = "item"
        self.item = None
        self.list = input_list

    def update_items(self):
        self.item = QGraphicsRectItem(
            0, 0, numbered_box_width * len(self.list), numbered_box_height
        )  # parent item: bkg
        # self.item.setBrush(QBrush(QColor(self.gene.color())))
        for boxdata_index in range(len(self.list)):
            boxdata = self.list[boxdata_index]
            if boxdata:
                number, hexcolor, filter = boxdata
                if filter == "filtered": # Should be only filtered
                    if number !=0:
                        colored_box = QGraphicsRectItem(
                            numbered_box_width * boxdata_index,
                            0,
                            numbered_box_width,
                            numbered_box_height,
                        )
                        colored_box.setBrush(QBrush(QColor(hexcolor)))
                        text_in_box = QGraphicsSimpleTextItem()
                        # printerr("***"+str(number), 1)
                        text_in_box.setText(str(number))
                        # Adding line 1
                        cross_brick_1 = QGraphicsLineItem(numbered_box_width* boxdata_index, numbered_box_height, numbered_box_width*boxdata_index+numbered_box_width,0)
                        pen = QPen(QColor("#F40000"))
                        pen.setWidth(3)
                        cross_brick_1.setPen(pen)
                        # Adding line 2
                        cross_brick_2 = QGraphicsLineItem(numbered_box_width*boxdata_index, 0, numbered_box_width*boxdata_index+numbered_box_width, numbered_box_height)
                        pen = QPen(QColor("#F40000"))
                        pen.setWidth(3)
                        cross_brick_2.setPen(pen)

                        font = text_in_box.font()
                        font.setPointSize(opt["nsize"])
                        text_in_box.setFont(font)  # setting default text size
                        reduce_font_if_necessary(
                            text_in_box,
                            numbered_box_width - 1,
                            numbered_box_height,
                            category="numbered box",
                        )
                        text_in_box.setZValue(
                            1
                        )  # default is 0, this is to be sure it is on top of colored box
                        text_in_box.setPos(boxdata_index * numbered_box_width + 1, 0)
                        colored_box.setParentItem(self.item)
                        text_in_box.setParentItem(self.item)
                        cross_brick_1.setParentItem(self.item)
                        cross_brick_2.setParentItem(self.item)
                # unfiltered
                else:
                  if number!=0:
                    colored_box = QGraphicsRectItem(
                        numbered_box_width * boxdata_index,
                        0,
                        numbered_box_width,
                        numbered_box_height,
                    )
                    colored_box.setBrush(QBrush(QColor(hexcolor)))
                    text_in_box = QGraphicsSimpleTextItem()
                    # printerr("***"+str(number), 1)
                    text_in_box.setText(str(number))
                    font = text_in_box.font()
                    font.setPointSize(opt["nsize"])
                    text_in_box.setFont(font)  # setting default text size
                    reduce_font_if_necessary(
                        text_in_box,
                        numbered_box_width - 1,
                        numbered_box_height,
                        category="numbered box",
                    )
                    text_in_box.setZValue(
                        1
                    )  # default is 0, this is to be sure it is on top of colored box
                    text_in_box.setPos(boxdata_index * numbered_box_width + 1, 0)
                    colored_box.setParentItem(self.item)
                    text_in_box.setParentItem(self.item)


    def _width(self):
        return self.item.rect().width()

    def _height(self):
        return self.item.rect().height()


class GeneFace(faces.Face):
    """ Colored rectangle for a gene, the color represent the label, the width and position how it spans the profile, the black or white lines the introns"""

    def __init__(self, limited_p2ghit_instance):
        faces.Face.__init__(self)
        self.type = "item"
        self.item = None
        self.gene = limited_p2ghit_instance

    def update_items(self):
        offset_for_additional_here = offset_for_additional * len(
            self.gene.additional_features
        )
        self.item = QGraphicsRectItem(
            0,
            0,
            gene_brick_width + offset_for_id + offset_for_additional_here,
            gene_brick_height,
        )  # parent item: bkg
        self.item.setPen(QPen(QtCore.Qt.NoPen))  # no black border
        if not opt["no_id"]:
            ### rectangle for prediction id
            chrom_id_rect = QGraphicsRectItem(
                0, gene_brick_height / 6, offset_for_id, gene_brick_height * 2 / 3
            )  # parent item: bkg
            chrom_id_rect.setParentItem(self.item)
            ## text for prediction id
            obj = QGraphicsSimpleTextItem()
            obj.setText(self.gene.id)
            reduce_font_if_necessary(
                obj, offset_for_id, gene_brick_height * 2 / 3, category="prediction id"
            )
            obj.setPos(
                (offset_for_id - obj.boundingRect().width()) / 2,
                (gene_brick_height - obj.boundingRect().height()) / 2,
            )  # centering the text
            obj.setParentItem(chrom_id_rect)

        for index, x in enumerate(self.gene.additional_features):
            # write(self.gene.full_id+' drawing secis '+' '+str(offset_for_additional), 1)
            secis_rect = QGraphicsRectItem(
                gene_brick_width + offset_for_id + offset_for_additional * index,
                gene_brick_height / 8,
                offset_for_additional - 1,
                gene_brick_height * 6 / 8,
            )  # parent item: bkg
            secis_rect.setParentItem(self.item)
            secis_rect.setBrush(QBrush(QColor(x.color())))
            # pen=QPen(); pen.width=3
            # pen.setColor(QColor("#FFFFFF"))
            # secis_rect.setPen(QPen(color="#FFFFFF", width=2) )
            obj = QGraphicsSimpleTextItem()
            obj.setText(x.text)
            obj.setBrush(QBrush(QColor("#FFFFFF")))
            reduce_font_if_necessary(
                obj,
                offset_for_additional_here,
                gene_brick_height * 6 / 8,
                category="secis description",
            )
            obj.setParentItem(secis_rect)
            obj.setPos(
                gene_brick_width
                + offset_for_id
                + offset_for_additional * index
                + (offset_for_additional - obj.boundingRect().width()) / 2,
                (gene_brick_height - obj.boundingRect().height()) / 2,
            )  # centering the text

        # drawing line to represent the 100 coverage for profile
        line_for_full_coverage = QGraphicsLineItem(
            offset_for_id, 0, offset_for_id + gene_brick_width, 0
        )
        line_for_full_coverage.setParentItem(self.item)
        ## rectangle for gene brick
        gene_brick = QGraphicsRectItem(
            offset_for_id + self.gene.relative_boundaries()[0] * gene_brick_width,
            0,
            (self.gene.relative_boundaries()[1] - self.gene.relative_boundaries()[0])
            * gene_brick_width,
            gene_brick_height,
        )
        gene_brick.setBrush(QBrush(QColor(self.gene.color())))
        gene_brick.setParentItem(self.item)

        if self.gene.label == "Missing":
            gene_brick = QGraphicsRectItem(offset_for_id, 0, gene_brick_width, gene_brick_height)
            gene_brick.setBrush(QBrush(QColor(self.gene.color())))
            gene_brick.setParentItem(self.item)

        # drawing lines for introns
        if not opt["I"] and len(self.gene.exons) > 1:
            cds_so_far = 0
            tot_cds = self.gene.length()
            for exon_index in range(len(self.gene.exons[:-1])):
                start, end = self.gene.exons[exon_index]
                cds_so_far += end - start + 1
                aa_so_far = cds_so_far / 3.0
                if self.gene.strand == "+":
                    intron_length = (
                        self.gene.exons[exon_index + 1][0] - end - 1
                    )  #######
                elif self.gene.strand == "-":
                    intron_length = (
                        start - self.gene.exons[exon_index + 1][1] - 1
                    )  #######
                if not self.gene.label == "Missing":
                    x_intron = (
                        self.gene.relative_position_in_ali_of(aa_so_far) * gene_brick_width
                        + offset_for_id
                    )
                    line_intron = QGraphicsLineItem(
                        x_intron, 1, x_intron, gene_brick_height - 2
                    )

                    color = opt["intron_color"]  # '#FFFFFF' #white
                    if intron_length <= 5:
                        color = "#EE0000"  # red for frameshifts
                    line_intron.setPen(QPen(QColor(color)))
                    line_intron.setZValue(1)
                    line_intron.setParentItem(gene_brick)

        # Adding cross for filtered species
        if opt["filter"]:
          if self.gene.filter == "filtered":
            # Adding line 1
            cross_brick_1 = QGraphicsLineItem(offset_for_id, 0, offset_for_id+gene_brick_width, gene_brick_height)
            pen = QPen(QColor("#F40000"))
            pen.setWidth(5)
            cross_brick_1.setPen(pen)
            cross_brick_1.setParentItem(gene_brick)
            # Adding line 2
            cross_brick_2 = QGraphicsLineItem(offset_for_id, gene_brick_height, offset_for_id+gene_brick_width, 0)
            pen = QPen(QColor("#F40000"))
            pen.setWidth(5)
            cross_brick_2.setPen(pen)
            cross_brick_2.setParentItem(gene_brick)



        if opt["f"]:
            for feature_name in self.gene.graphical_features:
                for aa_position in self.gene.graphical_features[feature_name]:
                    x_feature = (
                        self.gene.relative_position_in_ali_of(aa_position)
                        * gene_brick_width
                        + offset_for_id
                    )
                    line_feature = QGraphicsLineItem(
                        x_feature, 1, x_feature, gene_brick_height - 2
                    )
                    if feature_name in graphical_features_colors:
                        line_feature.setPen(
                            QPen(QColor(graphical_features_colors[feature_name]))
                        )
                    line_feature.setZValue(1.1)
                    line_feature.setParentItem(gene_brick)

        if not opt["T"]:
          if self.gene.filter != "filtered":
            ## text for chromosome
            obj = QGraphicsSimpleTextItem()
            obj.setPos(offset_for_id + 1, 0)
            obj.setText(self.gene.chromosome)
            font = obj.font()
            font.setPointSize(opt["bsize"])
            obj.setFont(font)  # setting default text size
            reduce_font_if_necessary(
                obj, gene_brick_width, category="text for chromosome"
            )
            obj.setZValue(2)
            obj.setParentItem(gene_brick)
            ## text for positions
            obj = QGraphicsSimpleTextItem()
            obj.setPos(offset_for_id + 1, gene_brick_height / 2 + 1)
            obj.setText(
                join([str(i) for i in self.gene.boundaries()], self.gene.strand)
            )  # e.g. 1+101 (positive strand) 45-80 (negative strand
            font = obj.font()
            font.setPointSize(opt["bsize"])
            obj.setFont(font)  # setting default text size
            reduce_font_if_necessary(obj, gene_brick_width, category="positions text")
            obj.setZValue(2.1)
            obj.setParentItem(gene_brick)

    def _width(self):
        return self.item.rect().width()

    def _height(self):
        return self.item.rect().height()


def is_selenoprofiles_title(title):
    """ Returns True if it is a selenoprofiles2 title, False if not """
    if (
        "chromosome:" in title
        and "strand:" in title
        and "positions:" in title
        and title.split()[0].count(".") in [6,5,4,2]
    ):
        return True
    return False


class gene_attribute(object):
    """ """

    def __init__(self, text=None, type=None, color=None):
        self.type = type
        self.text = text
        self.data = {}
        self.custom_color = color

    def color(self):
        if self.custom_color:
            return self.custom_color
        else:
            return "#000000"


class limited_p2ghit(gene):
    """This class is analog to p2ghit, but lacks some of its data. """

    def load_from_header(self, header):
        gene.load_from_header(self, header)
        if "species:" in header:
            tline = header.split("species:")[1]
            if tline.startswith('"'):
                species_name = tline[1:].split('"')[0]
            else:
                species_name = tline.split()[0]
        elif len(self.id.split(".")) >= 5:
            species_name = unmask_characters(
                replace_chars(self.id.split(".")[3], "_", " ")
            )
        else:
            species_name = "None"
        self.species = species(species_name)
        # self.program= header.split('prediction_program:')[1].split()[0]
        self.label = self.id.split(".")[2]
        self.filter = self.id.split(".")[5]
        self.profile_name = self.id.split(".")[0]
        if len(self.id.split(".")) >= 5:
            self.target_name = self.id.split(".")[-1]
        else:
            self.target_name = base_filename(self.target)
        self.full_id = self.id
        self.id = self.id.split(".")[1]
        self.relative_boundaries_data = (
            []
        )  # filled when the relative_boundaries function is called
        self.graphical_features = {}
        for feature_name in graphical_features_names:
            if feature_name + ":" in header:
                tt = header.split(feature_name + ":")[1]
                if tt and not tt[0] in " \n":
                    self.graphical_features[feature_name] = []
                    for block in tt.split()[0].split(","):
                        if "-" in block:
                            n, leng = list(map(int, block.split("-")))
                        else:
                            n = int(block)
                            leng = 1
                        self.graphical_features[feature_name].append(b, leng)
                    # float(n)
        self.additional_features = []

    def summary(self):
        return gene.summary(self, other_fields=["program", "label", "profile_name"])

    def relative_boundaries(self):
        """Return a list of two float numbers, indicating the boundaries in percentage of where the prediction is spanning the profile. max: 0.000001, 1.0 """
        if not self.relative_boundaries_data:
            try:
                self.profile
                assert self.profile
            except:
                raise Exception(
                    "relative_boundaries ERROR the .profile attribute must be defined to call this method! it failed on this object: "
                    + str(self)
                )
            self.relative_boundaries_data = [
                self.relative_position_in_ali_of(1),
                self.relative_position_in_ali_of(len(nogap(self.seq))),
            ]
        return self.relative_boundaries_data

    def relative_position_in_ali_of(self, position):
        """ Returns the position in the alignment corresponding to position pos in the aminacid sequence of self. Analog to position_in_ali of class alignment, but specific for this class. input is 1 based, output is a float, max=1.0 ... like relative boundaries"""
        pos_seq = 0
        for p, aa in enumerate(self.seq):
            if aa != "-":
                pos_seq += 1
            if pos_seq >= position:
                return float(p + 1) / len(self.seq)

    def color(self):
        """Return the hex representation of the color with which this gene will be drawn, depending on the label of this object. """
        if self.label in label_to_color:
            return label_to_color[self.label]
        else:
            return label_to_color["unknown"]

    def target_name(self):
        """Return the filename of the target for this prediction, removing the extension """
        return join(base_filename(self.target).split(".")[:-1], ".")

    def sequence_identity_with_profile(self, dont_count_gaps=0):
        """This maps the prediction back to the profile and returns the average sequence identity with the sequences of the profile. dont_count_gaps is flag for the following behavior:                                                                                                                                                                             
    0 -> terminal gaps are not counted                                                                                                                                           
    1 -> gaps are not counted                                                                                                                                                    
    2 -> everything is counted (resulting score will be much lower than the other methods for uncomplete hits)                                                                   
  """
        dont_count_gaps = int(dont_count_gaps)
        dcg = bool(dont_count_gaps == 1)
        dctg = bool(dont_count_gaps == 0)
        return self.alignment_with_profile(
            dont_shrink=True, title="UnMaTcHaBlE"
        ).average_sequence_identity_of(
            title="UnMaTcHaBlE", dont_count_gaps=dcg, dont_count_terminal_gaps=dctg
        )

    def alignment_with_profile(self, profile_ali="", dont_shrink=False, title=""):
        """ This functions is designed to be equivalent to the one with the same name in the p2ghit class of selenoprofiles, but to be working for limited_p2ghit
    """
        if not title:
            title = self.output_id()
        a = self.profile.copy()
        a.add(title, self.seq)
        a.remove_empty_columns()
        return a

    def output_id(self):
        return (
            self.profile_name
            + "."
            + self.id
            + "."
            + self.label
            + "."
            + self.species.name
            + "."
            + self.target_name
        )


graphical_features_names = []
graphical_features_colors = {}
gene_brick_width = 150
gene_brick_height = 40
offset_for_id = 25
offset_for_additional = 25

opt = {}


def main(args):
    #########################################################
    ############ loading options
    # global opt; opt=command_line(def_opt, help_msg, '*', synonyms=command_line_synonyms )
    # global temp_folder; temp_folder=Folder(random_folder(opt['temp'])); test_writeable_folder(temp_folder, 'temp_folder'); set_MMlib_var('temp_folder', temp_folder)
    global opt
    opt = args

    write("", 1)

    if opt['colors']:
        write('\nOption -colors: displaying the built-in colors for selenoprofiles labels\n', 1)
        for a, b, c in sorted(color_tuples, key=lambda x:x[0]):
            write(f'{a:<20}: {c:<16} {b}', 1)
        write('\nExiting...', 1)
        sys.exit()
    
    if not opt['i']:
        raise notracebackException("selenoprofiles drawer ERROR you must provide at least one .ali file with -i\nRun with -h to see usage")
    
    if opt["img_h"] == -1:
        opt["img_h"] = None

    if opt["img_w"] == -1:
        opt["img_w"] = None

    global gene_brick_width
    gene_brick_width = opt["W"]
    global gene_brick_height
    gene_brick_height = opt["H"]
    global numbered_box_width
    numbered_box_width = 30
    global numbered_box_height
    numbered_box_height = 40

    global offset_for_id
    if opt["no_id"]:
        offset_for_id = 0
    global offset_for_additional

    # checking input
    global tree_input_file
    if opt["t"]:
        tree_input_file = opt["t"]
        write(f"Reading input file: {tree_input_file}", 1)
        check_file_presence(tree_input_file, "tree_input_file")
        t = PhyloTree(tree_input_file)
        if opt["m"]:
            for node in t.traverse():
                node.name = unmask_characters(replace(node.name, "_", " "))
        for node in t:
            if node.is_leaf():
                node.is_used = 0
                node.columns = {}  # indexed with family name
    else:
        print(
            "WARNING no tree was provided in input: an unstructured tree with all species encountered in the input alignment will be used. To derive the phylogenetic tree of species included in ncbi taxonomy, visit: https://github.com/jhcepas/ncbi_taxonomy "
        )
        t = PhyloTree()

    tree_style = TreeStyle()
    if opt["C"]:
        write("Circular tree mode on", 1)
        tree_style.mode = "c"
        # tree_style.scale *= 10
        tree_style.allow_face_overlap = True
    else:
        tree_style.mode = "r"
    tree_style.branch_vertical_margin = 12
    tree_style.draw_aligned_faces_as_table = True
    tree_style.aligned_table_style = 1
    tree_style.show_leaf_name = False

    global explicit_labels
    explicit_labels = opt["explicit_labels"].split(",")
    global graphical_features_names
    global graphical_features_colors
    if opt["f"]:
        for piece in opt["f"].split(","):
            splt = piece.split("#")
            gf_name = splt[0]
            if len(splt) > 1:
                graphical_features_colors[gf_name] = "#" + splt[1]
            graphical_features_names.append(gf_name)

    write("", 1)
    global p2ghit_by_name
    p2ghit_by_name = {}
    global families_order
    families_order = []
    global labels_seen_for_family_hash
    labels_seen_for_family_hash = {}
    for ali_file in opt["i"]:
        if ":" in ali_file and "-" in ali_file.split(":")[1]:
            ali_file_position_range = [
                int(i) for i in ali_file.split(":")[1].split("-")
            ]  # [start, stop]
            ali_file = ali_file.split(":")[0]

        print("Input alignment: " + ali_file)
        ali = alignment(ali_file)

        profile_ali = alignment()
        for title in ali.titles():
            if not is_selenoprofiles_title(title):
                profile_ali.add(title, ali.seq_of(title))

        if not profile_ali.nseq():
            profile_ali.add("puppet", "X" * profile_ali.length())
        # profile_ali.remove_useless_gaps()

        any_title_was_included = False
        for title in ali.titles():
            # print "Title: "+title
            seq = ali.seq_of(title)
            seq_no_terminal_gaps = seq
            while seq_no_terminal_gaps.startswith("-"):
                seq_no_terminal_gaps = seq_no_terminal_gaps[1:]
            while seq_no_terminal_gaps[-1] == "-":
                seq_no_terminal_gaps = seq_no_terminal_gaps[:-1]

            if is_selenoprofiles_title(title):
                any_title_was_included = True
                x = limited_p2ghit()
                x.load_from_header(title)
                x.seq = seq
                x.profile = profile_ali
                family = x.profile_name
                profile_ali.name = family  ## will be the same for each selenoprofiles title in this alignment . nonetheless we don't want to know it before here otherwise we'd need to make assumptions

                species_name = x.species.name

                p2ghit_by_name[x.full_id] = x
                if family not in labels_seen_for_family_hash:
                    labels_seen_for_family_hash[family] = {}
                if x.label not in labels_seen_for_family_hash[family]:
                    labels_seen_for_family_hash[family][x.label] = 1

                    if opt["a"] and opt["c"] and not x.label in explicit_labels:
                        labels_seen_for_family_hash[family]["others"] = 1

                if not opt[
                    "t"
                ]:  # no tree was specified. Let's build a puppet tree with the organisms names
                    try:
                        node = t & species_name
                    except:
                        node = t.add_child(name=species_name)
                        node.is_used = 0
                        node.columns = {}  # indexed with family name
                        node = t & species_name
                try:
                    node = [
                        n for n in t.search_nodes(name=species_name) if n.is_leaf()
                    ][0]
                    node.is_used = 1

                    if not node.is_leaf():
                        raise Exception(
                            "ERROR species: "
                            + species_name
                            + " is not a leaf in the input tree!"
                        )
                    if family not in node.columns:
                        node.columns[family] = []
                    node.columns[family].append(x)

                except IndexError:
                    #          print species_name
                    #          raise
                    if not opt["g"]:
                        raise Exception(
                            "ERROR can't find a node in the tree for species: "
                            + species_name
                        )
                    else:
                        print("ignoring species not found: " + species_name)

        if any_title_was_included:
            sorted_l = sorted(list(labels_seen_for_family_hash.keys()))
            for family in sorted_l:
              if family not in families_order:
                families_order.append(family)
                title_face = faces.TextFace(family, fsize=opt["tsize"])
                title_face.hz_align = 1
                tree_style.aligned_header.add_face(title_face, column=len(families_order))
        else:
            printerr(
                "WARNING no single title was included for alignment: " + ali_file, 1
            )

    write("", 1)
    if opt["add"]:
        write(f'Reading additional gene features from {opt["add"]}', 1)
        for line in open(opt["add"]):
            stripped = line.strip()
            if stripped:
                splt = stripped.split("\t")
                p2gname = splt[0]
                text = replace(splt[1], "\\n", "\n")
                if len(splt) > 2 and splt[2]:
                    color = splt[2]
                else:
                    color = None

                p2ghit_by_name[p2gname].additional_features.append(
                    gene_attribute(text=text, color=color)
                )

    if opt["sp_add"]:
        write(f'Reading species features from {opt["sp_add"]}', 1)
        for line in open(opt["sp_add"]):
            stripped = line.strip()
            if stripped:
                splt = stripped.split("\t")
                species = splt[0].strip()
                try:
                    node = t & species
                except:
                    try:
                        species = unmask_characters(replace(species, "_", " "))
                        node = t & species
                    except:
                        raise Exception(
                            "ERROR can't find a node in the tree for species in -sp_add file: "
                            + str([species])
                        )

                    value = splt[1]
                    if value.startswith("COLOR"):
                        node.species_coloring = "#" + value.split("#")[1]
                    else:
                        if not hasattr(node, "species_attributes"):
                            node.species_attributes = []
                        node.species_attributes.append(value)

    # parse tree and prune useless nodes!
    nodes_to_keep = [t]
    for node in t:
        if node.is_used:
            nodes_to_keep.append(node)

    n_to_prune = len([1 for n in t]) - len(nodes_to_keep) + 1
    if n_to_prune:
        if not opt["e"]:
            write(f"Pruning {n_to_prune} species with no associated results", 1)
            t.prune(nodes_to_keep)
        else:
            write(
                "Option -e: keeping {n_to_prune} species with no associated results", 1
            )

    # setting outgroup
    if opt["o"]:
        if opt["o"] == "":  # in [1, True]:
            write("Setting outgroup automatically with get_midpoint_outgroup", 1)
            t.set_outgroup(t.get_midpoint_outgroup())
        elif opt["o"] == "None":
            write("Skipping setting outgroup", 1)
        else:
            write(f'Setting outgroup: {opt["o"]}', 1)
            outgroup = str(opt["o"])
            matches = t.search_nodes(name=outgroup)
            # print "*********"+str(matches)
            if len(matches) > 0:
                t.set_outgroup(matches[0].name)
            else:
                raise Exception(
                    "ERROR the species "
                    + outgroup
                    + " was not found in the tree."
                    + " Maybe it is because this species had no results and was pruned out: try with option -e"
                    * int(not opt["e"])
                )

    write("Ladderizing tree", 1)
    t.ladderize()
    # print t
    write("", 1)
    if not opt["out"]:
        write("Opening interactive tree visualization ...", 1)
        t.show(mylayout, tree_style)
    else:
        write(f'Writing to file --> {opt["out"]}', 1)
        t.render(
            opt["out"],
            layout=mylayout,
            tree_style=tree_style,
            h=opt["img_h"],
            w=opt["img_w"],
        )

    for category in warnings_reduced_text:
        printerr(
            "WARNING "
            + str(warnings_reduced_text[category])
            + " "
            + category
            + " text(s) were reduced in point size to fit the dedicated space",
            1,
        )

    if opt["prompt"]:
        interactive_mode(message="Tree is loaded into variable: t")

    ###############


def mylayout(node):

    ## Phylogeny layout
    if node.is_leaf():
        ## Setting the leaf color name

        # node.img_style['bgcolor']='#DDDDDD'
        column_index = 0

        ## modify to have a customized species name printed
        main_species = node.name
        main_species_in_filenames = replace_chars(
            mask_characters(main_species), " ", "_"
        )
        if opt["common"]:
            main_species = common_names.setdefault(
                main_species_in_filenames, main_species_in_filenames
            )

        if hasattr(node, "species_coloring"):
            fgcolor = node.species_coloring
        else:
            fgcolor = "#000000"

        nameFace = faces.TextFace(
            main_species, fgcolor=fgcolor, fsize=opt["ssize"]
        )  # species name!
        nameFace.vt_align = 1
        faces.add_face_to_node(
            nameFace, node, column=column_index, position="branch-right"
        )
        column_index += 1

        if hasattr(node, "species_attributes"):
            for value in node.species_attributes:
                valueFace = faces.TextFace(value, fgcolor="#000000", fsize=opt["ssize"])
                valueFace.vt_align = 1
                faces.add_face_to_node(
                    valueFace, node, column=column_index, position="branch-right"
                )
            column_index += 1

        for family in families_order:

            if opt["a"]:
                ####################### ABSTRACT MODE
                if family in node.columns:
                    if opt["F"]:
                        familyNameFace = faces.TextFace(
                            family, fsize=opt["bsize"], fgcolor="#000000"
                        )
                        familyNameFace.margin_left = 20
                        familyNameFace.margin_right = 20
                        faces.add_face_to_node(
                            familyNameFace, node, column=column_index, aligned=True
                        )

                    count_per_label = {}
                    for gene_index in range(len(node.columns[family])):
                        x = node.columns[family][gene_index]
                        key = (x.label, x.filter)
                        if key not in count_per_label:
                            count_per_label[key] = 0
                        count_per_label[key]+=1

                    # create rectangle with width of: box_width * len(labels_seen_for_family_hash[family].keys())
                    list_to_build_numbered_box = []
                    if not opt["c"]:
                        for label in ordered_labels:
                            if label in labels_seen_for_family_hash[family]:
                                key_unfiltered = (label, "unfiltered")
                                key_filtered = (label, "filtered")
                                if key_unfiltered not in count_per_label and key_filtered not in count_per_label:
                                    list_to_build_numbered_box.append([])
                                else:
                                    count_unfiltered = count_per_label.get(key_unfiltered, 0)
                                    count_filtered = count_per_label.get(key_filtered, 0)
                                    list_to_build_numbered_box.append(
                                        [
                                            count_unfiltered,
                                            label_to_color.setdefault(
                                                label, label_to_color["unknown"]
                                            ),
                                            "unfiltered",
                                        ]
                                    )
                                    list_to_build_numbered_box.append(
                                        [
                                            count_filtered,
                                            label_to_color.setdefault(
                                                label, label_to_color["unknown"]
                                            ),
                                            "filtered",
                                        ]
                                    )
                    else:  # condensating boxes in max n columns. possible labels are defined by explicit_labels option

                        count_per_label["others"] = 0
                        for label in ordered_labels:
                            if (
                                label in count_per_label
                                and not label in explicit_labels
                            ):
                                count_per_label["others"] += count_per_label[label]
                        for label in explicit_labels + ["others"]:
                            if label in labels_seen_for_family_hash[family]:
                                if label not in count_per_label:
                                    list_to_build_numbered_box.append([])
                                elif not count_per_label[label]:
                                    list_to_build_numbered_box.append([])
                                else:
                                    list_to_build_numbered_box.append(
                                        [
                                            count_per_label[label],
                                            label_to_color.setdefault(
                                                label, label_to_color["unknown"]
                                            ),
                                        ]
                                    )
                    numbered_boxes_face = NumberedBoxFace(list_to_build_numbered_box)
                    numbered_boxes_face.margin_right = opt["margin_boxes"]
                    numbered_boxes_face.hz_align = 1
                    numbered_boxes_face.rotable = False
                    faces.add_face_to_node(
                        numbered_boxes_face, node, column=column_index, aligned=True
                    )

                else:
                    if opt["F"]:
                        emptyFamilyNameFace = faces.TextFace(family, fgcolor="#FFFFFF")
                        emptyFamilyNameFace.margin_left = 4
                        emptyFamilyNameFace.margin_right = 2
                        faces.add_face_to_node(
                            emptyFamilyNameFace, node, column=column_index, aligned=True
                        )

            else:
                ####################### NORMAL MODE
                if family in node.columns:
                    if opt["F"]:
                          familyNameFace = faces.TextFace(family, fgcolor="#000000")
                          familyNameFace.margin_left = 2
                          familyNameFace.margin_right = 2
                          faces.add_face_to_node(
                              familyNameFace, node, column=column_index, aligned=True
                          )

                    list_of_genes = node.columns[family]
                    try:
                          list_of_genes.sort(
                              key=lambda x: ordered_labels.index(x.label)
                          )  # this list will host the ordered list of genes to draw in this species. The ordered is determined by the appearances of labels in ordered_labels
                    except:
                          for g in list_of_genes:
                              if not g.label in ordered_labels:
                                  print(
                                      " WARNING label "
                                      + g.label
                                      + " was not found among the known ones. "
                                  )

                    for gene_index in range(len(list_of_genes)):
                          x = list_of_genes[gene_index]
                          gene_face = GeneFace(x)
                          gene_face.margin_left = 5
                          gene_face.margin_right = 5
                          faces.add_face_to_node(
                              gene_face, node, column=column_index, aligned=True
                          )

                    #            if opt['ali']:

                    # add separator?

                    else:
                      if opt["F"]:
                          emptyFamilyNameFace = faces.TextFace(family, fgcolor="#FFFFFF")
                          emptyFamilyNameFace.margin_left = 2
                          emptyFamilyNameFace.margin_right = 2
                          faces.add_face_to_node(
                              emptyFamilyNameFace, node, column=column_index, aligned=True
                          )

            column_index += 1


#######################################################################################################################################
def prune_tree(t, nodes_to_keep):
    """ Fixed (and faster) prunning algorithm. Use this until I fix
    the problem within the main ETE branch.

    'nodes_to_keep' must be the list of node instances that you want
    to keep in the final tree. All nodes must be leaves, if not, they
    are automatically converted into leaves by removing their
    children.

    So far, this function is quite verbose. Printing slows down a bit
    the process, but you can follow the progress...
    """
    # print "Getting tree path..."
    # Converts to set to speed up searches
    if type(nodes_to_keep) == set:
        to_keep = nodes_to_keep
    else:
        to_keep = set(nodes_to_keep)

    # print "Checking that all nodes are leaves..."
    not_leaves = [n for n in nodes_to_keep if not n.is_leaf()]
    if len(not_leaves) > 0:
        #    print "\nFixing", len(not_leaves), "non-leaf nodes..."
        # Converts all internal species nodes into leaves by removing all
        # their sub-species or strains
        for nl in not_leaves:
            for c in nl.get_children():
                c.detach()
    to_detach = []

    # print "Selecting unused nodes"
    counter = 0
    for node in t.traverse("postorder"):
        #    print "\r", counter,
        counter += 1
        for c in node.children:
            if c in to_keep:
                to_keep.add(node)
                break
        if node not in to_keep:
            to_detach.append(node)
            for c in node.children:
                to_detach.remove(c)
    # print "\nDetaching", len(to_detach), "nodes"
    counter = 0
    for node in to_detach:
        #    print "\r", counter,
        counter += 1
        node.detach()
    # print "\nFixing", len(to_keep), "orphan nodes"
    counter = 0
    for node in to_keep:
        #    print "\r", counter,
        counter += 1
        if len(node.children) == 1:
            node.delete()

    if len(t.children) == 1:
        try:
            a = t.children[0]
            a.delete()

        except:
            pass

    ############################

    return t


##########################3


# def close_program():
#   if opt['debug']: raw_input('check temp folder:'+temp_folder)
#   if 'temp_folder' in globals() and is_directory(temp_folder):
#     bbash('rm -r '+temp_folder)
#   try:
#     if get_MMlib_var('printed_rchar'):
#       printerr('\r'+printed_rchar*' ' ) #flushing service msg space
#   except:
#     pass


# if __name__ == "__main__":
#   try:
#     main()
#     close_program()
#   except Exception:
#     close_program()
#     raise

# else:
#   global opt;
#   opt=get_MMlib_var('opt')
