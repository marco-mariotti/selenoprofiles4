from string import *
import sys
from subprocess import *

#sys.path.insert(0, "/home/mmariotti/software/selenoprofiles/libraries/")
#sys.path.append("/home/mmariotti/software/selenoprofiles")

from .MMlib3 import *

# from profiles_classes import *

help_msg = """selenoprofiles join: collects results obtained on different targets and joins them into a single file per profile.

This program looks for .ali files, which are normally produced by selenoprofiles (i.e. its -output_ali option is on by default)
You may directly provide the list of .ali files (Usage #2), or let the program automatically find them (Usage #1).
In all cases, the source profile.fa files are loaded from the -profiles_folder

Usage #1:  selenoprofiles join -d selenoprofiles_output_folder  -o selenoprofiles_join_output 
Usage #2:  selenoprofiles join -i list_of_ali_files             -o selenoprofiles_join_output 

* Options (the + indicates an argument is required):
-o  +   output folder of this program, where it will create joined .ali files

* Options for Usage #1:
-d  +  folder previously used as output of selenoprofiles, containing results for multiple targets
-p  +  profile(s) to be searched. For accepted arguments, see selenoprofiles -h. If not provided, searches all builtin profiles
-pf +  file with list of profile(s), one per line. Overrides -p

* Options for Usage #2:
-i  +   file with list of input .ali files

* Other options
-da       do not alter selenoprofiles ids. Normally the species and target name is added to distinguish between results on different targets
-u        do no transfer alignments: the output is unaligned with predictions and profile sequences. Increases speed
-ds       do not shrink the alignment. Normally, the program detects gappy columns due to the alignment transfer procedure, and realigns them with mafft
-config   path to selenoprofiles config file, used to load default options
-h OR --help    print this help and exit """

# command_line_synonyms = {
#     "p": "fam",
#     "p_list": "fam_list",
#     "P": "fam",
#     "profile": "fam",
#     "s": "species",
#     "s_list": "species_list",
#     "t": "target",
#     "t_list": "target_list",
# }

def_opt = {  #'temp':'/home/mmariotti/temp',
    "i": 0,
    "u": 0,
    "ds": 0,
    "o": "",
    "profile": "",  #has syn p
    'pf':"",
    "d": "",
    'da':0,
    "debug": 0,
}


#########################################################
###### start main program function


def is_selenoprofiles_output_title(title):
    return (
        "chromosome:" in title
        and "target:" in title
        and "positions:" in title
        and "strand:" in title
    )


def main(output_folder,
         fam2filelist=None,
         opt={}):
    #########################################################
    ############ loading options
    # global opt
    # if not args:
    #     opt = command_line(def_opt, help_msg, "do", synonyms=command_line_synonyms)
    # else:
    #     opt = args
    #set_MMlib_var("opt", opt)
    # global temp_folder; temp_folder=Folder(random_folder(opt['temp'])); test_writeable_folder(temp_folder, 'temp_folder'); set_MMlib_var('temp_folder', temp_folder)
    # global split_folder;    split_folder=Folder(opt['temp']);               test_writeable_folder(split_folder); set_MMlib_var('split_folder', split_folder)
    
    ################### file list (attempts) is now defined. loading alignments and transferring them to get a joined ali
    for fam in sorted(fam2filelist.keys()):
        joined_ali = False
        for file_attempt in fam2filelist[fam]:

            # if opt["debug"]:
            #     write("ATTEMPT: " + file_attempt, 1)
            if is_file(file_attempt):
                write("Reading: " + file_attempt, 1)
                this_ali = alignment(file_attempt)
                corrected_ali = alignment()
                for title in this_ali.titles():
                    if is_selenoprofiles_output_title(title):

                        g = gene()
                        g.load_from_header(title)
                              
                        if g.species:
                            new_title = (
                                g.id
                                + "."
                                + replace_chars(mask_characters(g.species), " ", "_")
                                + "."
                                + join(base_filename(g.target).split(".")[:-1], "_")
                            )
                        else:
                            raise Exception("ERROR missing species argument in title: "+title+" from file "+file_attempt)
                        dont_print_until_double_commas = False
                        for i in title.split()[1:]:
                            if dont_print_until_double_commas:
                                dont_print_until_double_commas = i[-1] != '"'
                            else:
                                if not i.startswith("species:"):
                                    new_title += " " + i
                                elif i.startswith('species:"'):
                                    dont_print_until_double_commas = True

                    else:
                        new_title = title
                    corrected_ali.add(new_title, this_ali.seq_of(title))

                if opt["u"]:
                    if not joined_ali:
                        joined_ali = alignment()
                    for new_title in corrected_ali.titles():
                        if not joined_ali.has_title(new_title):
                            joined_ali.add(
                                new_title, nogap(corrected_ali.seq_of(new_title))
                            )

                else:  #### normal case: now transfering the alignment
                    if joined_ali:
                        joined_ali = joined_ali.transfer_alignment(
                            corrected_ali, dont_shrink=True
                        )
                    else:
                        joined_ali = corrected_ali

        if joined_ali:
            outfile = output_folder + fam + ".ali"
            # if opt["suffix"]:
            #     outfile = output_folder + fam + "." + opt["suffix"] + ".ali"
            if not opt["u"] and not opt["ds"]:
                joined_ali.shrink()
            write("Writing ", how='magenta')
            write("---------> " + outfile, 1)
            joined_ali.display(outfile)


#######################################################################################################################################


# def close_program():
#     if "temp_folder" in globals() and is_directory(temp_folder):
#         bbash("rm -r " + temp_folder)
#     try:
#         if get_MMlib_var("printed_rchar"):
#             printerr("\r" + printed_rchar * " ")  # flushing service msg space
#     except:
#         pass


# if __name__ == "__main__":
#     try:
#         main()
#         close_program()
#     except Exception:
#         close_program()
#         raise
