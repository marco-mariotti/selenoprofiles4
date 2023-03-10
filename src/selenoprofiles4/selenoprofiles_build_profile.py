#!/home/mmariotti/miniconda3/envs/selenoprofiles/bin/python -u
from string import *
import sys
from subprocess import *
import obonet

from .MMlib3 import *
from .selenoprofiles4 import (
    profile_alignment,
    get_selenoprofiles_var,
    load,
    parse_blast,
)


help_msg = (
    """selenoprofiles build : build an alignment into a profile that can be used to search genomes. Usage:  

selenoprofiles build -i input_alignment.fa -o built_profile.fa   [-attribute "value"] [options]

The basics:
- Only aligned fasta (.fa or .fasta) is accepted as input. Preferably, use an alignment of 10 sequences or more
- The default options will work to create a suitable configuration for most use cases
- The input filename determines the profile name, which can only contain alphanumeric characters or underscores

* I/O options:
-i          input file in aligned fasta file format.  If -o is not provided, it is overwritten with output
-o          output file with extension .fa or .fasta. The filename determines the profile name. Further output files are created using this prefix
-load ARG   load a .config file (ARG) from an existing profile. Any options on command line will overwrite the configuration loaded from this file
-y          avoid interactive prompting for confirmation. Only available if -o is provided

* Pre-processing alignment:
-r  [ARG]   remove sequence redundancy before building the profile, leaving a max grade of identity between any sequence pair which is ARG [default 0.9]
-t  [ARG]   before building, the alignment is trimmed removing the columns with more than ARG gaps [default 0.1 or number of sequences minus one]

* Inspecting profile features:
-p          open the interactive environment after building (or loading) the profile
-d  [ARG]   draw the distribution of profile AWSIc scores with pylab interactive graphics. Use 2 as ARG to plot AWSIw, 3 as ARG to plot average sequence identities
-D  [ARG]   print the conservation score for each sequence. Arguments are like for -d
-out  ARG   drawn output file if option -d is active. Accepted formats are pdf and png

* Profile configuration:
"""
    + profile_alignment.make_profile.__doc__
    + """

* System options
-temp  +          temporary folder. A folder with random name is created here, used and deleted at the end of the computation.
-print_opt        print currently active options
-h OR --help      print this help and exit
"""
)


'''Pending development:
* GO tools:
-go         prints GO terms found in input sequences. It runs a blast search against Uniref50 
-GO         set the GO terms in the profile. This is not recommended unless you check them first

'''

command_line_synonyms = {}

def_opt = {
    "temp": "",
    "load": 0,
    "i": 0,
    "o": 0,
    "v": 0,
    "Q": 0,
    "c": 0,
    "o": "",
    "r": "",
    "t": "",
    "GO": 0,
    "go": 0,
    "C": 0,
    "p": 0,
    "d": 0,
    "out": 0,
    "y":0,
}


#########################################################
###### start main program function


def main(opt, selenoprofiles_config):
    #########################################################
    ############ loading options
    #global opt
    # if not args:
    #     opt = command_line(
    #         def_opt, help_msg, "io", synonyms=command_line_synonyms, nowarning=1
    #     )
    # else:
    #opt = args
    set_MMlib_var("opt", opt)
    # config_filename = opt["sp_config"]

    global temp_folder
    temp_folder = Folder(random_folder(opt["temp"]))
    test_writeable_folder(temp_folder, "temp_folder")
    set_MMlib_var("temp_folder", temp_folder)
    # global split_folder;    split_folder=Folder(opt['temp']);               test_writeable_folder(split_folder); set_MMlib_var('split_folder', split_folder)
    # checking input
    global input_file
    input_file = opt["i"]
    check_file_presence(input_file, "input_file")
    
    if opt["d"]:
        try:
            import pylab
        except:
            raise Exception(
                "ERROR import pylab failed: pylab modules must be installed to use option -d\nInstall it by running:  conda install matplotlib"
            )
    keep_awsi = False
    if opt["d"] or opt["D"]:
        if opt["d"] and opt["D"] and opt["d"] != opt["D"]:
            raise Exception(
                "ERROR both option -d and -D are specified but they are different! exiting..."
            )
        if opt["d"]:
            conservation_score_flag = opt["d"]
        elif opt["D"]:
            conservation_score_flag = opt["D"]
        if conservation_score_flag in [1, 2]:
            keep_awsi = True

    ## loading selenoprofiles configuration file. the only thing we're interested in is the gi2go_db
    # selenoprofiles_config = configuration_file(config_filename)

    if opt["go"] or opt["GO"]:
        # try:
        #     import annotations.GO.Parsers.oboparser as oboparser
        # except:
        #     raise Exception, "ERROR the gene ontology utilities must be installed to use option -go or -GO! Please see selenoprofiles installation script."
        uniref2go_db = selenoprofiles_config["uniref2go_db.DEFAULT"]
        uniref_db = selenoprofiles_config["tag_db.DEFAULT"]
        GO_obo_file = selenoprofiles_config["GO_obo_file"]
    ####

    # preparing build_options that will be passed to make_profile method of profile_alignment
    build_options = opt.copy()
    # if build_options.has_key('queries'):  exec( "build_options['queries']="+build_options['queries'])  #coverting queries to appropriate type (list of integers or list of strings)
    if (
        "go_terms" in build_options
    ):  # coverting go_terms to appropriate type (list of strings)
        ss = build_options["go_terms"]
        build_options["go_terms"] = []
        if ss[0] == "[":
            ss = ss[1:]
        if ss[-1] == "]":
            ss = ss[:-1]
        for s in ss.split(","):
            if s[0] in "'\"":
                s = s[1:]
            if s[-1] in "'\"":
                s = s[:-1]
            build_options["go_terms"].append(s)

    if opt["load"]:
        # loading old config file
        check_file_presence(
            opt["load"], "profile configuration provided with option -load"
        )
        if opt["load"].split(".")[-1] != "config":
            raise Exception(
                "ERROR the file provided with option -load should have the extension .config"
            )
        old_profile_configuration = configuration_file(opt["load"])
        for k in old_profile_configuration:
            if k not in build_options:
                build_options[k] = old_profile_configuration[k]

    if opt['y'] and not opt['o']:
        raise Exception("ERROR option -y is accepted only if output is specified with -o")
    
    if not opt['y']:
        write('\nReady to build profile with configuration:\n', 1)
        m='\n'.join([f'{k} = {build_options[k]}'  for k in build_options  if not k in def_opt and not k.startswith('__')])
        if not m:
            write('<all profile configuration is default>', 1)
        else:
            write(m, 1)

        z='\n'.join([f'{k} = {build_options[k]}'  for k in build_options  if k in def_opt and build_options[k]!=def_opt[k] and  not k in ['temp', 'i', 'o']])
        if z:
            write('\nMoreover, this/these procedure(s) will be carried out:', 1)
            write(z, 1)
                    
        write('', 1)

        write( 'The profile will be written to files:')
        ff=input_file if not opt['o'] else opt['o']
        write(f'\n  --> {ff}'  + (' WARNING overwriting an existing file!' if is_file(ff) else ''  ))
        write(f'\n  --> {ff}.config'  + (' WARNING overwriting an existing file!' if is_file(ff+'.config') else ''  ))
        write(f'\n  --> {ff}.profile_data'  + (' WARNING overwriting an existing file!' if is_file(ff+'.profile_data') else ''  ))
        
        confirm=None
        while confirm!='Y':
            confirm=input('\n\nEnter Y to confirm, or Enter N or Ctrl+C to cancel\n>')
            if confirm=='N':
                write('Cancelled, exiting...', 1)
                sys.exit()    

    if opt["o"]:
        # preparing output file, if specified
        print("Copying " + input_file + " to " + opt["o"] + " ...")
        bbash("cp " + input_file + " " + opt["o"])
        input_file = opt["o"]

    if opt["r"]:
        # removing redundancy
        if opt["r"] == 1:
            opt[
                "r"
            ] = 0.9  # setting default value for max pairwise identity after trimming
        print(
            "Removing redundancy of "
            + input_file
            + " to max pairwise identity "
            + str(opt["r"])
            + "... "
        )
        a = alignment(input_file)
        a.remove_redundancy(opt["r"], inplace=True, silent=False)
        a.display(input_file)

    if opt["t"]:
        a = alignment(input_file)
        # trimming columns
        if opt["t"] == 1:
            opt["t"] = 0.1
            if a.nseq() < 1 / opt["t"]:
                opt["t"] = 1 / float(a.nseq()) - 0.01
        print("Trimming columns with threshold: " + str(opt["t"]) + "... ", end=" ")
        old_length = a.length()
        a.trim_columns(max_non_gaps=opt["t"], inplace=True, remove_empty_seqs=True)
        new_length = a.length()
        print(str(old_length - new_length) + " columns removed")
        a.display(input_file)

    if opt["GO"] or opt["go"]:
        # determining the GO codes annotated for the proteins in the profile

        # loading gene ontology
        #parserO = oboparser.Parser(open(GO_obo_file))
        #gene_ontology = parserO.parse()

        print('Loading gene ontology ... ')
        gene_ontology = obonet.read_obo(GO_obo_file) #2023
        #gene_ontology.define_annotations_level(silent=1)

        representative_seqs = temp_folder + "representative_seqs.fa"
        a = alignment(input_file)
        repre_titles = a.titles()
        random.shuffle(repre_titles)
        repre_titles = repre_titles[:5]
        write_to_file(
            join([">" + t + "\n" + nogap(a.seq_of(t)) for t in repre_titles], "\n"),
            representative_seqs,
        )

        blast_out = input_file + ".blast_uniref"
        temp_blast = temp_folder + "blast_uniref"
        blast_cmnd = 'blastall -i {i} -o {o} -d {d} -b 50 -v 50 -p blastp -F "m S" -f 999  -e 1e-10 -M BLOSUM80 -G 9 -E 2; mv {o} {f}'.format(
            i=representative_seqs, o=temp_blast, d=uniref_db, f=blast_out
        )
        if not is_file(blast_out):
            print(
                "Running blastp of 5 random sequences in your alignment against Uniref50 ..."
            )
            bbash(blast_cmnd)

        # putting in a file the id for the profile proteins
        id_list_file = temp_folder + "id_list_file"
        id_codes_list = [bh.chromosome.split()[0] for bh in parse_blast(blast_out)]
        write_to_file(join(id_codes_list, "\n"), id_list_file)

        # determining gi_GO associations
        print("Obtaining GO of the sequences matched by blast ...")
        id_go_associations_hash = {}
        id_go_associations_string = bbash(
            "gawk -v id_file="
            + id_list_file
            + """ -F"\\t" 'BEGIN{ while ((getline idline < id_file)>0){ GI_INPUT[idline]=1 } } { split($1, GI, "; "); split("", gi_match); for (i=1; i<=length(GI); i++){ if (GI[i] in GI_INPUT) { gi_match[GI[i]]=1}   }; o=""; for (g in gi_match) o=o"; " g; if (o)  print substr(o, 3) "\\t" $2  }' """
            + uniref2go_db,
            dont_die=1,
        )

        if id_go_associations_string:
            for line in id_go_associations_string.split("\n"):
                for id_code in line.split("\t")[0].split(
                    "; "
                ):  # more than one gi can be present in a line
                    id_go_associations_hash[id_code] = []
                    for go_term in line.split("\t")[1].split("; "):
                        # keeping only those for molecular function
                        # if 'GO:0003674' in gene_ontology.get_all_parents_ids(go_term):      ## done by construction of the id2go file
                        id_go_associations_hash[id_code].append(go_term)

        # counting how many times each GO is found
        go_count_hash = {}
        for go_list in list(id_go_associations_hash.values()):
            for go_code in go_list:
                if go_code not in go_count_hash:
                    go_count_hash[go_code] = 0
                go_count_hash[go_code] += 1

        interactive_mode('see gene_ontology and go_count_hash')

        # adding count for the children
        for go_code in go_count_hash:
            try:
                go_term_obj = gene_ontology.get_term_by_id(go_code)
                child_go_terms_list = gene_ontology._internal_dag.predecessors(go_code)
                for child_go_term in child_go_terms_list:
                    if child_go_term in go_count_hash:
                        go_count_hash[go_code] += go_count_hash[child_go_term]
            except IOError: # annotations.GO.ontology.NoSuchTermError:
                printerr(
                    "GO utility WARNING GO code not found: "
                    + go_code
                    + " ; this may be obsolete and it was removed from the GO database in "
                    + GO_obo_file
                    + " , or the database may need update. "
                )

        print("GO id     | N in ali | N in all prots | GOname")
        for go_code in go_count_hash:
            # cycling all gos found in profile
            go_term_obj = gene_ontology.get_term_by_id(go_code)
            child_go_terms_list = gene_ontology._internal_dag.predecessors(
                go_code
            )  # quite hidden a method in the ontology class, I had to browse to get it. hope it is right
            go_to_count_file = temp_folder + "go_to_count_file"
            write_to_file(join([go_code] + child_go_terms_list, "\n"), go_to_count_file)

            count_this_go_and_children = int(
                bbash(
                    "grep -c -w -f " + go_to_count_file + " " + uniref2go_db, dont_die=1
                )
            )
            print(
                go_term_obj.id.ljust(15)
                + str(go_count_hash[go_code]).ljust(10)
                + str(count_this_go_and_children).ljust(14)
                + go_term_obj.name
            )

        # removing those which have a parent considered
        all_go_codes_considered = list(go_count_hash.keys())
        for go_code in all_go_codes_considered:
            all_parents_of_this_go = gene_ontology.get_all_parents_ids(go_code)
            for parent_go in all_parents_of_this_go:
                if parent_go in all_go_codes_considered:
                    print(
                        "Removing "
                        + go_code
                        + " since its parent "
                        + parent_go
                        + " is already in the list of considered GOs."
                    )
                    del go_count_hash[go_code]
                    break

        if opt["go"]:
            print(
                'GO check  *** the proposed GO terms are:    -go_terms ["'
                + join(list(go_count_hash.keys()), '","')
                + '"]'
            )
            print("now exiting (without building the profile)...")
            sys.exit()

        if opt["GO"]:
            print(
                'Setting GO terms for this profile: ["'
                + join(list(go_count_hash.keys()), '","')
                + '"]'
            )
            build_options["go_terms"] = list(go_count_hash.keys())

    build_options["keep_awsi"] = keep_awsi

    ########### building!
    #  print build_options
    try:
        load(
            config_filename, {"puppet_option": 1, "temp": temp_folder}
        )  # putting a fake option just to give a non-null hash as args; this will avoid selenoprofiles trying to attempt loading options form commandline
    except:
        pass  # the load function crashes since it doesn't find the file necessary for run, for example the target file. nonetheless the variables are in place.

    print("Building profile ...")
    a = profile_alignment(input_file, build=True, **build_options)  # building profile!
    print(a.summary())

    print("*** " + a.filename + ".config successfully built! Here's its content:\n")
    print(bbash("cat " + a.filename + ".config"))

    if opt["c"]:
        print(
            "\nThis is the list of the fasta titles labelled as queries:\n"
            + join(a.queries_titles(), "\n")
        )

    if opt["d"] or opt["D"]:
        if opt["d"] and opt["D"] and opt["d"] != opt["D"]:
            raise Exception(
                "ERROR both option -d and -D are specified but they are different! exiting..."
            )
        if opt["d"]:
            conservation_score_flag = opt["d"]
        elif opt["D"]:
            conservation_score_flag = opt["D"]

        if conservation_score_flag == 3:
            title2score = {}  ## filling it now
        elif conservation_score_flag == 1:
            title2score = a.conservation_data["awsi_scores_with_coverage"]
        elif conservation_score_flag == 2:
            title2score = a.conservation_data["awsi_scores_without_coverage"]

        for title in a.titles():
            if title == "BLAST_QUERY_MASTER":
                continue
            if conservation_score_flag == 3:
                title2score[title] = a.average_sequence_identity_of(title)
            if opt["D"]:
                write(
                    title.split()[0].ljust(80)
                    + " "
                    + str(round(title2score[title], 3)),
                    1,
                )

        ordered_titles = list(title2score.keys())
        ordered_titles.sort(key=title2score.get)
        scores = []
        for title in ordered_titles:
            scores.append(title2score[title])

        average_score = average(scores)
        std_dev = std_deviation(scores)

        if opt["d"]:
            n_bins = 20
            n, bins, patches = pylab.hist(
                scores, bins=n_bins, normed=0, facecolor="green", alpha=0.75
            )
            y = max(n) * 0.6
            y2 = y * 0.9
            y3 = y * 0.6
            y4 = y * 0.2
            pylab.plot([average_score], [y], "D", color="red")
            pylab.plot(
                [average_score - std_dev, average_score + std_dev],
                [y2, y2],
                "D",
                color="purple",
            )
            pylab.plot(
                [average_score - 2 * std_dev, average_score + 2 * std_dev],
                [y3, y3],
                "D",
                color="blue",
            )
            pylab.plot(
                [average_score - 3 * std_dev, average_score + 3 * std_dev],
                [y4, y4],
                "D",
                color="aqua",
            )

            if opt["d"] == 1:
                title = "AWSIc scores of " + base_filename(input_file)
            if opt["d"] == 2:
                title = "AWSIw scores of " + base_filename(input_file)
            if opt["d"] == 3:
                title = "Average seq identities of " + base_filename(input_file)
            pylab.title(title)
            #      pylab.xlabel('Seq identity')
            pylab.ylabel("Frequency count")
            pylab.grid(True)
            pylab.xlim(0.0, 1.0)
            pylab.ylim(min(n), max(n))
            if opt["out"]:
                print(
                    "Writing image with score distribution to file: " + str(opt["out"])
                )
                pylab.savefig(opt["out"])  ## saving output file
            else:
                print("Opening pylab interactive graphics with score distribution ... ")

                pylab.show()  ## opening pylab intereactve environment

    if opt["p"]:
        interactive_mode(message="Profile alignment is loaded into variable: a")


#######################################################################################################################################


def close_program():
    if "temp_folder" in globals() and is_directory(temp_folder):
        bbash("rm -r " + temp_folder)
    try:
        if get_MMlib_var("printed_rchar"):
            printerr("\r" + printed_rchar * " ")  # flushing service msg space
    except:
        pass


if __name__ == "__main__":
    try:
        main()
        close_program()
    except Exception:
        close_program()
        raise
