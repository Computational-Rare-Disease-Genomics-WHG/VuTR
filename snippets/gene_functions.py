# gene_functions.py
# J. Martin

# pylint: skip-file
# flake8: noqa

from gene_classes import Gene


def reset_vars(vars_dict):
    vars_dict["Temp_Transcript_ID"] = ''
    vars_dict["Temp_Chr"] = ''
    vars_dict["Temp_Gene_Start_Coord"] = ''
    vars_dict["Temp_Gene_End_Coord"] = ''
    vars_dict["Temp_Gene_Strand"] = ''
    vars_dict["Temp_Gene_ID"] = ''
    vars_dict["Temp_Gene_Symbol"] = ''
    vars_dict["Temp_HGNC_Symbol"] = ''
    vars_dict["Temp_NCBI_Symbol"] = ''
    vars_dict["Temp_Gene_Description"] = ''
    vars_dict["Temp_cDNA_Sequence"] = ''
    return vars_dict


def read_MANE(filename):
    FirstEntry = True

    # Temporary variables
    tempvars = {
        "Temp_Transcript_ID": '',
        "Temp_Chr": '',
        "Temp_Gene_Start_Coord": '',
        "Temp_Gene_End_Coord": '',
        "Temp_Gene_Strand": '',
        "Temp_Gene_ID": '',
        "Temp_Gene_Symbol": '',
        "Temp_HGNC_Symbol": '',
        "Temp_NCBI_Symbol": '',
        "Temp_Gene_Description": '',
        "Temp_cDNA_Sequence": '',
    }

    Genelist = []
    # Stores all Gene objects as the file is read

    with open(filename, 'r') as f:
        # for-else clause executed when the loop terminates through exhaustion of the list
        gene_count = 0
        for line in f:
            if line[0] == '>':  # ie, new gene reached

                if FirstEntry == True:
                    # If this is the first time a gene is read, then no previous information collected, Parse header and proceed

                    gene_parse(line, tempvars)
                    FirstEntry = False  # Only reached after first gene recognised, this in future iterations will pass through if statement 2

                else:  # elif(FirstEntry == False):  # Note, elif here?
                    # All data, including sequence for previous data stored in temp variables, thus instantiate object
                    # with this data, reset, and begin to store information for new gene
                    gene_count += 1
                    print(gene_count)
                    print(tempvars["Temp_Gene_Symbol"])
                    tempvars["Temp_Gene_Symbol"] = Gene(tempvars)
                    # Construct Object, initialising with temporary variables
                    Genelist.append(tempvars["Temp_Gene_Symbol"])
                    # - Add object to data structure

                    reset_vars(tempvars)  # - Reset Variables
                    gene_parse(
                        line, tempvars
                    )  # - Parse Header and assign new variables

            if line[0] != '>':
                tempvars["Temp_cDNA_Sequence"] += line.rstrip()

        else:
            gene_count += 1
            print(gene_count)
            print(tempvars["Temp_Gene_Symbol"])
            tempvars["Temp_Gene_Symbol"] = Gene(tempvars)
            Genelist.append(tempvars["Temp_Gene_Symbol"])
            reset_vars(tempvars)

    return Genelist


def gene_parse(line, tempvars):

    # Add formal definition so can be found in "help"
    """
    The following code essentially reads in the data in a header file, parses it using
    layers of split functions, and assigns the correct info to temporary variables, which will
    be used to instantiate an Gene object
    """

    Header_Data = line.split(" ", 7)
    # Header Data - Splits Header into 8 chunks, where the splits occur at ' '

    Pre_Parse_Transcript_ID = Header_Data[0].split(">")
    # Pre_Parse_Transript_ID - Splits first chunk (eg >ENST00000342066.8) into two, sepeperated by '>', and takes the 2nd element
    tempvars["Temp_Transcript_ID"] = Pre_Parse_Transcript_ID[1]

    Pre_Parse_Chr_Coord_Strand = Header_Data[2].split(":")
    # Pre_Parse_Chr_Coord_Strand - Further splits third "Header Data" chunk (chromosome:GRCh38:1:925731:944574:1 gene:ENSG00000187634.12)
    # into smaller parts, seperated by ":" - these parts are then assigned to the relevant information
    tempvars["Temp_Chr"] = Pre_Parse_Chr_Coord_Strand[2]
    tempvars["Temp_Gene_Start_Coord"] = Pre_Parse_Chr_Coord_Strand[3]
    tempvars["Temp_Gene_End_Coord"] = Pre_Parse_Chr_Coord_Strand[4]
    tempvars["Temp_Gene_Strand"] = Pre_Parse_Chr_Coord_Strand[5]

    Pre_Parse_Gene_ID = Header_Data[3].split(":")
    # Pre_Parse_Gene_ID - Further splits fourth "Header Data" chunk (eg gene:ENSG00000187634.12) into two,
    # seperated by ":", where the second element is used for the gene ID
    tempvars["Temp_Gene_ID"] = Pre_Parse_Gene_ID[1]

    Pre_Parse_Gene_Symbol = Header_Data[6].split(":")
    # Pre_Parse_Gene_Symbol = splits seventh chunk (eg gene_symbol:SAMD11), across ":", and so symbol can be assigned
    tempvars["Temp_Gene_Symbol"] = Pre_Parse_Gene_Symbol[1]
    # if Pre_Parse_Gene_Symbol[1] == "SMIM40":
    #     print ("TEST NUMBER")
    #     sys.exit()

    Pre_Parse_Gene_Description, Pre_Parse_Gene_Info = Header_Data[7].split("[", 1)
    # Pre_Parse_Gene_Description = splits eighth chunk (eg description:sterile alpha motif domain containing 11 [Source:HGNC Symbol;Acc:HGNC:28706])
    # into 2, seperated by "[", takes the first element, further splits by ":", and takes the second element, before cutting off whitespace

    # Pre_Parse_Gene_Info = splits eighth chunk (eg description:sterile alpha motif domain containing 11 [Source:HGNC Symbol;Acc:HGNC:28706])
    # into 2, seperated by "[", takes the second element, further splits according to ":", takes the last element, and removes the end "]"

    Pre_Parse_Symbol_Check = Pre_Parse_Gene_Info.split(" ", 1)[0]
    # Pre_Parse_Symbol_Check - splits up chunk in order to see the source of the Gene information as HGNC / NCBI

    if Pre_Parse_Symbol_Check == "Source:HGNC":
        # If the source is HGNC, this block is entered, the test is split, and temporary variables are assigned
        Pre_Parse_HGNC_Symbol = Pre_Parse_Gene_Info.split(":")[3]
        tempvars["Temp_HGNC_Symbol"] = (
            "HGNC:" + Pre_Parse_HGNC_Symbol[:-2]
        )  # remove end "]"

    elif Pre_Parse_Symbol_Check == "Source:NCBI":
        # If the source is NCBI, this block is entered, the test is split, and temporary variables are assigned
        Pre_Parse_NCBI_Symbol = Pre_Parse_Gene_Info.split(":")[2]
        tempvars["Temp_NCBI_Symbol"] = Pre_Parse_NCBI_Symbol[:-2]  # remove end "]"

    Pre_Parse_Gene_Description = Pre_Parse_Gene_Description.split(":")[1]
    tempvars["Temp_Gene_Description"] = Pre_Parse_Gene_Description[
        :-1
    ]  # remove end whitespace

    return tempvars


def compile_gene_structure(filename, gene_search):

    with open(filename, 'r') as f:
        unsorted_gene_structure = (
            []
        )  # Build the initial list which will store gene segments
        gene_name = ""  # Initialise gene_name variable
        transcript_name = ""
        gene_found = False  # Turns true when match is found, which will break out of program after structure is compiled

        for line in f:
            if (
                line[0] == 'c'
            ):  # Only reads appropriate lines containing gene information
                line_info = line.split("\t")  # Spits the values seperated by tabs

                if line_info[2] == "gene":
                    if (
                        gene_found == True
                    ):  # If this condition is met, the intended gene has been found previously, and further search is unnecessary

                        segment_count = 0
                        loop_count = 0
                        if (
                            strand == "+"
                        ):  # If Forward strand, sorting will be different

                            # Defines new list containing the segments compiled below, but sorted in ascending co-ordinate order
                            gene_structure = sorted(
                                unsorted_gene_structure, key=lambda k: k['Beginning']
                            )

                            # Adds Start and End indexes for each segment, so that each sequence index mapped to genetic co-ordinates
                            for elem in gene_structure:
                                if segment_count == 0:
                                    indexes = {
                                        "Beginning_Index": 0,
                                        "End_Index": (elem["Length"] - 1),
                                    }
                                    elem.update(indexes)

                                elif segment_count != 0:
                                    start_ind = (gene_structure[segment_count - 1])[
                                        "End_Index"
                                    ] + 1
                                    indexes = {
                                        "Beginning_Index": start_ind,
                                        "End_Index": (start_ind + elem["Length"] - 1),
                                    }
                                    elem.update(indexes)

                                # print(elem)
                                segment_count += 1
                                loop_count += 1

                            return gene_structure

                        elif strand == "-":
                            # Defines new list containing the segments compiled below, but sorted in descending co-ordinate order (due to reverse strand direction)
                            gene_structure = sorted(
                                unsorted_gene_structure,
                                key=lambda k: k['Beginning'],
                                reverse=True,
                            )

                            # Adds Start and End indexes for each segment, so that each sequence index mapped to genetic co-ordinates
                            for elem in gene_structure:
                                if segment_count == 0:
                                    indexes = {
                                        "Beginning_Index": 0,
                                        "End_Index": (elem["Length"] - 1),
                                    }
                                    elem.update(indexes)

                                elif segment_count != 0:
                                    start_ind = (gene_structure[segment_count - 1])[
                                        "End_Index"
                                    ] + 1
                                    indexes = {
                                        "Beginning_Index": start_ind,
                                        "End_Index": (start_ind + elem["Length"] - 1),
                                    }
                                    elem.update(indexes)

                                # print(elem)
                                segment_count += 1
                                loop_count += 1
                            return gene_structure

                    # Parse line to extract Gene ID and Symbol
                    pre_gene_info = line_info[8].split(";")
                    gene_name = ((pre_gene_info[3].split("="))[1])[:-1]
                    gene_id = (pre_gene_info[0].split("="))[1]
                    strand = line_info[6]
                    # print(gene_name)
                    # print(gene_id)
                    # print(strand)

                # Parse line to extract Transcript ID
                if line_info[2] == "transcript":
                    pre_gene_info = line_info[8].split(";")
                    transcript_name = (pre_gene_info[0].split("="))[1]
                    # print(transcript_name)

                if (
                    (gene_search == gene_name)
                    or (gene_search == gene_id)
                    or (gene_search == transcript_name)
                ):
                    # If item matches the Gene name, gene ID or transcript ID searched for my the user
                    gene_found = True

                    # Conditional entry based on lines that describe a gene segment
                    if (
                        (line_info[2] == "CDS")
                        or (line_info[2] == "five_prime_UTR")
                        or (line_info[2] == "three_prime_UTR")
                    ):
                        if (
                            strand == "+"
                        ):  # Fw/Rv strand determines assignment of beginning and end
                            # Defines Dictionary values based on the parsed information in the file
                            gene_segment = {
                                "Type": line_info[2],
                                "Beginning": line_info[3],
                                "End": line_info[4],
                                "Length": (int(line_info[4]) - int(line_info[3]) + 1),
                            }
                            # Defines a list with containing segments, which are to be ordered
                            unsorted_gene_structure.append(gene_segment)

                        elif strand == "-":
                            # Defines Dictionary values based on the parsed information in the file
                            gene_segment = {
                                "Type": line_info[2],
                                "Beginning": line_info[4],
                                "End": line_info[3],
                                "Length": (int(line_info[4]) - int(line_info[3]) + 1),
                            }
                            # Defines a list with containing segments, which are to be ordered
                            unsorted_gene_structure.append(gene_segment)


def reading_frame(seq, init_index):
    # Read 3 Bases from initial position, shifting 3 bases per iteration - index values
    # range from initial index, to end of sequence-2, in steps of 3
    return [
        (seq[index : index + 3]) for index in range(init_index, len(seq) - 2, 3)
    ]  # Return series of base triplets


# Read 3 Bases from initial position, shifting 3 bases per iteration - index values
# range from initial index, to end of sequence-2, in steps of 3


def generate_frames(seq):
    reading_frames = []
    reading_frames.append(reading_frame(seq, 0))
    reading_frames.append(reading_frame(seq, 1))
    reading_frames.append(reading_frame(seq, 2))
    return reading_frames


def calculate_Kozak_strength(sequence, start_index):
    strength = None
    # #  To annotate the strength of the Kozak consensus into which the uAUG was formed we assessed the positions at
    # #  −3 and +3 relative to the A of the AUG, known to be the most important bases for dictating strength of translation.
    # #  If both the −3 base was either A or G and the +3 was G, Kozak was annotated as ‘Strong’,
    # #  if either of these conditions was true, Kozak was deemed to be ‘Moderate’ and if neither was the case ‘Weak’
    plus3base = sequence[start_index + 3]
    minus3base = sequence[start_index - 3]

    if ((minus3base == "A") or (minus3base == "G")) and (plus3base == "G"):
        strength = "Strong"

    elif ((minus3base == "A") or (minus3base == "G")) or (plus3base == "G"):
        strength = "Moderate"

    else:
        strength = "Weak"

    return strength


def uORF_gene_region(gene_structure, index):

    region = None  # Placeholder for full title including co-ordinates
    region_name = None  # Placeholder for region
    for segment in gene_structure:
        if (index >= segment["Beginning_Index"]) and (index <= segment["End_Index"]):

            # Conditionals used to format output in more natural form
            if segment["Type"] == "five_prime_UTR":
                region_name = "5'UTR"
            elif segment["Type"] == "CDS":
                region_name = "CDS"
            elif segment["Type"] == "three_prime_UTR":
                region_name = "3'UTR"

            region = (
                region_name
                + " :: "
                + str(segment["Beginning"])
                + " - "
                + str(segment["End"])
            )
            return region


def generate_ORF_list(read_frames, n, sequence, gene_structure):
    # Must build in functionality to:
    # - Only calculate the uORFs where ATG is before the beginning of the coding sequence
    # - Calculate Region
    # - Calculate co-ordinates
    # - Calculate Kozak
    ORFs = []
    ORF = []
    codon_count = 0
    open_RF = False
    initial_start_found = False
    # open_RF condition set up so that stop codons will only be recognised if they are terminating a uORF

    for codon in read_frames[n]:
        # Iterate over each frame, storing each in frame start (+ relevant data) and in frame
        # stop (+ relevant data) in a seperate dictionary
        # The dictionaries are then added to a list uORF, which stores the initial start, and all in-frame starts until a stop is reached
        # These will then be added to a list which gives all uORFS in the gene

        if codon == "ATG":
            open_RF = True

            start_index = n + (
                codon_count * 3
            )  # Gets sequence index of the beginning of the start codon

            # print("START", start_index, (seq[start_index -3], seq[start_index +3]))
            # print(start_index)
            # print(seq[start_index], seq[start_index + 1], seq[start_index + 2])
            # print(seq[start_index -3], seq[start_index +3])

            kozak_strength = calculate_Kozak_strength(sequence, start_index)
            start_region = uORF_gene_region(gene_structure, start_index)
            # print(kozak_strength)
            print(
                "START",
                start_index,
                start_region,
                (sequence[start_index - 3], sequence[start_index + 3]),
                kozak_strength,
            )
            start = {  # Defines Dictionary values
                "type": "start",
                "index": start_index,
                "region": start_region,
                "kozak_strength": kozak_strength
                # <Coords>:<> - deal with this later
            }

            print(start)
            ORF.append(start)
            if (
                initial_start_found == False
            ):  # Captures index of ORF's first start, for later sorting
                orf_start = start_index
                orf_region = start_region
                initial_start_found = True

        if (open_RF == True) and (
            (codon == "TAA") or (codon == "TAG") or (codon == "TGA")
        ):
            # Condition only satisfied for correct codon existing with an ORF
            stop_index = n + (codon_count * 3)
            stop_region = uORF_gene_region(gene_structure, stop_index)
            print("STOP", stop_region, stop_index)

            open_RF = (
                False  # Resets the value to indicate that an ORF hasn't been started
            )
            stop = {"type": "stop", "index": stop_index, "region": stop_region}
            print(stop)

            ORF.append(stop)  # Add stop codon's dictionary to uORF list
            ORF.append(orf_region)
            ORF.append(orf_start)  # Adds ORF start to list, used when sorting
            ORFs.append(ORF)  # Add uORF list to the uORFs list
            ORF = []  # Reset uORF list for the next start codon
            initial_start_found = False
            orf_start = None
            orf_region = None

        # def get_coords():
        #     Get genetic co-ordinates of the
        codon_count += 1
    return ORFs
