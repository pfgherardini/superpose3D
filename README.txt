- Installation

In order to compile superpose from sources you need CMake. CMake is a free and
open source cross-platform build system available at http://www.cmake.org/.
Cmake will be used to create a build environment specific for the OS you are
using. Detailed instructions for each system are given in the following
paragraphs.



- Mac OS

In order to compile superpose (or any C/C++ program for that matter!) you need
Xcode. This is an optional install available in the DVD you have received when
you bought your computer. Xcode is also available at
http://developer.apple.com/technology/xcode.html


Decompress the superpose archive, go into the Applications folder and start
CMake. In the CMake gui point the "Where the source code is" field to the
location of the superpose source code folder. Similarly point the "Where to
build the binaries" to the build/ subdir in the superpose folder. Next press
the configure button (bottom part of the window) and select as generator "Unix
Makefiles".  Press configure again until there are no more values in red in
the main portion of the window. Afterwards press the generate button. Open the
terminal go into the build/ subdirectory and type
	
	make

The build products will be located in build/src. If you want to inspect/modify
the sources you can create an Xcode project by choosing Xcode as the
generator. The project will be created in the build/ subdirecoty. Double click
on the file and everything should open in Xcode.

Note: if you installed cmake via fink you can do all of the above from the
command line.  Decompress the superpose archive, go into the build/ subdir and
type

	cmake .. -G"Unix Makefiles" 
	make

The build products will be located in build/src. If you want to create an
Xcode project type
	
	cmake .. -GXcode

Important note: if you have troubles compiling the program, esepcially if you
get weird linking errors, you may have a version of Xcode that was not
installed in the version of the OS you are currently using. This happens for
instance when you migrate your programs (including Xcode) from an older Mac
with an earlier version of the OS. If this is the case please reinstall Xcode.


- Linux

You will need to have g++ installed. g++ will certainly be available in the
software repositories of your distribution so check there.
Decompress the superpose archive, go into the build/ subdir and type

	cmake .. 
	make

The build products will be located in build/src. 

- Windows

CMake can create projects for a variety of development environments in Windows
(e.g. Borland, Microsoft Visual C++, Cygwin, Eclipse etc.). The instructions
that follow refer to Visual C++ but the steps are similar with any development
environment. Visual C++ Express is available for free at
http://www.microsoft.com/express/download/ Decompress the superpose archive.
Navigate the windows "start" menu to CMake(Cmake-gui) and start the CMake gui.
In the CMake gui point the "Where the source code is" field to the location of
the superpose source code folder. Similarly point the "Where to build the
binaries" to the build/ subdir in the superpose folder. Next press the
configure button (bottom part of the window) and select as generator the
version of Visual C++ that you are using.  Press configure again until
there are no more values in red in the main portion of the window. Afterwards
press the generate button.  This will create a Visual C++ project inside the
build/ subdir. Double-click on the project to open it with Visual C++. In the
Visual C++ window set the build type to Release using the drop-down menu
located in the middle of the tool bar and navigate to build -> Build
ALL_BUILD. This will start the compiler. The build products will be located in
the Release subdirectory of build/src. If you have troubles with all this
please refer to the Visual Studio documentation.  BIG FAT WARNING: If you are
on Windows specify all the paths that you give to superpose using / as the
directory separator!

- Files included in the distribution

The following files are all located in the data/ subdir. The parameters file
we provide will allow you to start comparing structures right away, without
the trouble of developing a structure representation yourself.

pdb_res_dict.txt: a file containing all the aminoacids known to the program.
Please remember to set the pdb_res_dict option in the parameters files to
point to this file (see below, "Specifying options and parameters")

query3d_par.txt: a file of parameters corresponding to the representation used
in [Ausiello et al., BMC Bioinformatics. 2005 Dec 1;6 Suppl 4:S5].

pints_par.txt: a file of parameters corresponding to the representation used
in [Stark et al., Nucleic Acids Res. 2003 Jul 1;31(13):3341-4].

schmitt_par.txt: a file of parameters corresponding to the representation used
in [Schmitt et al., J Mol Biol. 2002 Oct 18;323(2):387-406]. 


- A bit of terminology

In the following documentation a "pseudoatom" is a point used by superpose to
represent a residue or part of it. A pseudoatom may correspond to an actual
protein atom, e.g. you could specifiy that a pseudoatom called CA is used to
represent the CA of a residue. However this is not always true. For instance
you could define a pseudoatom called "bar" which corresponds to the geometric
centroid of the side chain atoms. The term "pseudoatom" therefore designates
the points making up a residue in the internal representation used by the
program.  A "match" is a correspondence between the pseudoatoms describing two
protein structures. The match may involve complete residues, single
pseudoatoms, or a combination of the two (see below, "Modes"). Obviously,
since this is a pairing, a match will always involve the same number of
elements in the two structures at hand.

- Description of the algorithm

The structural comparison algorithm uses a branch & bound strategy to find the
largest subset of residues (or part of them) between two protein structures
that can be superimposed under a given RMSD threshold, irrespective of their
position along the sequence. The user is completely free to specify the way in
which protein structures will be represented during the comparison.  Users can
also specify equivalence criteria between residues. Only residues which have
been defined as equivalent will be paired.


The algorithm starts by creating all the possible pairs of equivalent
pseudoatoms in the two structures (i.e. matches of length 1). These matches
are extended by generating all the possible length 2 matches starting from
each length 1 match. For this purpose all the possible pairs of residues (or
pseduoatoms depending on the mode) which are close to those already belonging
to the match (i.e. having a pair of pseudoatoms at distance lower than a given
threhsold) are addedd to the first two. For instance, let the first aminoacid
in the probe structure have i neighbours and the corresponding matched
aminoacid in the target structure have j neighbours. The algorithm generates i
x j new matches of length 2. All the matches are evaluated and kept if their
constituent residues can be superimposed with a RMSD lower than the threshold
currently in use. The optimal superimposition between two sets of points is
calculated using the Quaternion method [Coutsias et al., J Comput Chem 2004,
25:1849-1857]. The matches that have been retained are then recursively
extended to length 3 and so on.

- Running time and complexity

The complexity of the procedure is clearly exponential in the size of the
probe and target chains.However using appropriate RMSD and given the average
size of protein chains, meaningful results can generally be obtained rather
quickly. The bottom line is that the algorithm runs faster as a greater
number of matches fail in the early stages of the exploration. This allows an
early pruning of the search tree which results in fewer matches having to be
considered for extension. If superpose is taking too long to complete you can speed it up
by lowering the RMSD threshold, and evenutally choosing stricter pairing
rules, so that fewer combinations of residues are considered. Also critical is the
structural similarity of the proteins given as input. The more similar they
are, the more seed matches that will have to be considered for extension. Please
remember that this software is about _local_ structural similarities. If you
have two proteins which are globally similar you will be better off using a
fold comparison method such as DALI [Holm and Sander, J Mol Biol 1993,
233:123-138].

- Command line syntax and files needed by superpose

You can launch the program using the following syntax

	superpose -i [input file] -p [parameters file] [-o [out file]] [-v]

The command line arguments can be specified in any order, -o is optional as
the name of the output file can also be provided as an option in the input
file (see below). The -v option turns on verbose output. In addition you will
also need the pdb_res_dict.txt file located in the data/ directory of the
superpose distribution. This is a file containing a list of all the residues
known to the program, along with their constituent atoms. If an aminoacid is
not in there it will not be considered when parsing the PDB file. The name of
these residues and their constituent atoms were taken from the Chemical
Component Dictionary (CCD) available from the PDB website
(http://www.wwpdb.org/ccd.html). The syntax of the file used by superpose,
which you are free to modify to add your custom residues if necessary, is as
follows:

	[3 letter code] [3 letter code of the parent amino acid] [one letter code of the parent aa] [","-separated list of heavy atoms]

The parent aminoacid is the standard aminoacid from which a modified aa
derives (e.g. for phospho-threonine, TPO, the parent aa is threonine, THR).
This corresponds to the _chem_comp.mon_nstd_parent_comp_id field in the CCD
file.


- Specifying options and parameters

Options and parameters used by superpose are always specified with the syntax

	[key] [value]

where key and value are separated by a single space. Each key/value pair has
to go on a different line except for the "def" and "equiv" keys that can span
multiple lines (see below). A description of the options and parameters
recognized by the program follows. The letter in parenthesis after the key
name specifies whether the option should be placed in the parameters (P) or
input (I) file. This letter is only used here and is not part of the key name.


	def (P): this is the keyword used to prefix lines containing the
	definition of protein residues (see "Specifying protein residues")
	
	equiv (P): this is the keyword used to prefix lines containing the
	definition of equivalences (see "Specifying equivalences")

	rmsd_threshold (P): the maxium RMSD of a match, default 0.7

	pdb_res_dict (P): the full path to the file containing the dictionary of
	PDB residues known to the program, see "Files needed by superpose"

	neighbour_threshold (P): the distance threshold under which two residues
	(or pseudoatoms depending on the mode, see "Modes") are considered
	neighbours, see "Description of the algorithm" for additional details.

	score_min (P): the minimum length of a match, default 3

	score_max (P): the maximum length of a match, default 10

	max_longest_matches (P): the maximum number of matches of maximum score
	that will be reported, after this limit is reached the program will
	stop looking for additional matches between the same pair of proteins.
	Please note: if you hit this limit the output you are looking at is
	very likely _not_ what you want. The matches reported by the program
	will be the _first_ found during the exploration. There might very
	well be better matches of the same length that are not reported
	because this limit was reached and therefore the exploration stopped. This
	threshold is essentially a safeguard against the possibility that
	superpose will waste a huge amount of time in the comparison of two
	closely related structures (see "Description of the algorithm" for
	clarifications). Default is 5; -1 means unlimited (Use at your own risk!).

	mode (P): the mode used by superpose (see, "Modes"). Default res; possible
	values: res, res_multiple, atom

	probe (I): the protein chain(s) to be used as probe(s). These are specified
	as follows, separated by ";":

		[pdb file]([chain][:residue range])

	residue ranges are specified either as a ","-separated list of residue
	numbers (as specified in the resSeq PDB field) or as a "-"-separated
	lower and upper range limits, or as a combination of both. The chain
	is mandatory. You can type * as the chain and the program will load
	all the chains contained in the file. Instead of a single file you can
	also specify the name of a directory by enclosing it in square
	brackets. In this case the program will traverse the directory
	recursively, trying to load all the files it finds. In this case, if
	you specify a chain or a residue range the expression will be applied
	to each file. Examples of valid probe definitions:

		probe pdb101m.ent(A:12,34,80-90,40-50);pdb102m.ent(A)
		probe [subdir/](*) #This means all the chains in all the files contained in the subdirectory "subdir"
	
	It is also possible to specify as input a zone comprising all the
	residues within a given distance from a central residue. The syntax is
	as follows:

		[pdb file]([chain][:central residue]z[distance threshold])

	For instance

		probe pdb101m.ent(A:120z6.5)

	The above statement selects all the residues whose distance from
	residue 120 of chain A is less than 6.5 Angstroms. To inspect which
	residues the program is selecting you can use the -v command line
	option for turning verbose output on.
	
	target (I): the protein chain(s) to be used as target(s). The syntax
	is the same as for probes. If only probes are provided the program
	will compare each one of them against all the others. Otherwise each
	probe will be compared with each target.

	pdb_dir (I): the directory containing the pdb_files, defaults to the
	current directory. If a directory is given this will be prepended to
	all the file and directory names specified in "probe" and "target"

	output_file(I): full path to the output file. Defaults to the basename
	of the input wit "_out.txt" appended. Warning: if the file already
	exists it will be overwritten!



- Specifying protein residues

Each residue is specified in a line starting with the keyword "def". The
format is as follows:
	
	def [PDB residue name]=[atom1]:[pseudoatom name];[atom2]:[pseudatom name]

PDB residue name should be the three-letter code of the residue as found in
the PDB file (e.g. ALA for alanine). The "=" sign is followed by a
";"-delimited list of atoms that define the residue at hand. The atom naming
scheme is the one defined in the PDB Chemical Compound Dictionary. In the case
of glycine residues the program will add a "fake" CB atom while parsing the PDB
file. The coordinates of this atom are calculated from the N, CA and C atoms
so as to have a tetrahedral geometry at the CA.  You can also specify a name
of your choice for each atom after the ":" (e.g. CA:bar, associates the name
bar to the CA atom). The rationale for assigning custom names will be clear
when we discuss residue equivalences. If you do not specify any custom name
the atom name will be used to label the point.


It is also possible to define a "pseudoatom" as the geometric centroid of a
user-specified list of PDB atoms.
The syntax is as follows:

	avg(CA,C,O):bar

This clause defines a pseudoatom named "bar" whose coordinates are calculated
as the geometric centroid of the CA, C and O atoms. Note that in this case it
is mandatory to assign a name to the pseudoatom (i.e. the above clause would
be invalid without the ":bar" part).  Inside the "avg" clause you can also use
the keyword "side_chain" which designates all the atoms except the main chain
ones (e.g. N, CA, C, O).

The definitions that follows should clarify the syntax.

	def ALA=CA;N;C
	def ALA=CA;N:bar;C

The above definitions all specify alanine as represented by three atoms (CA,
N, C). The second one also assigns the label "bar" to the N atom.

	def TRP=CA;avg(side_chain):foo

This defines TRP as represented by the CA plus a pseudoatom called "foo" and
calculated as the geometric centroid of the side chain.

	def GLU=avg(OE1,OE2,CG):acidic

This specifies GLU as the geometric centroid of the OE1, OE2, and CG atoms.
This pseudoatom is assigned the label "acidic"
Note that the same atom can appear in different avg clauses, e.g. the
following definition is perfectly valid

	def GLU=agv(OE1,CG):acidic1;avg(OE2,CG):acidic2

The only requirement is that the user-specified names must be unique in the
definition of each residue (e.g.: using "acidic" for both pseudoatoms would
have been invalid). However they may repeat between residues (e.g. two different
residues may each have a pseudoatom called "acidic").


- Specifying equivalences

Once you have defined the points that the program should use to represent
residues you have to provide a list of equivalences. These statements specify
which residues are allowed to match and which points should be used for the
superimposition. These deifinitions are given in lines starting with the
"equiv" keyword.
The syntax is as follows

	equiv [res]{.[atom]}=[res]{.[atom]}=[res]{.[atom]}

The .[atom] part is optional so let's start with the basic case.

	equiv ALA=GLY=VAL

Tells the program that the residues ALA, GLY and VAL are equivalent and can be
matched with each other. Note that if you want identical residues to match you
have to explicitly specify it! e.g. if the above line is the only equivalence
you provide the program will _not_ pair ALA with ALA, GLY with GLY and VAL with
VAL. If you want this three residues to pair with each other and also with
themselves you have to used the following lines

	equiv ALA=ALA
	equiv GLY=GLY
	equiv VAL=VAL
	equiv ALA=GLY=VAL

Even though this may seem redundant it allows you to exclude certain residue
pairs altogether. For instance if you do not want to see matches between
hydrophobic residues you can leave them out from the list of equivalences and
you will not even see matches of residues of the same type.

An alternative way to specify equivalences is the following

	equiv ALA=GLY;VAL

Note that the residues on the right side of the ‘=‘ sign are separeted by a ‘;’. 
This means that ALA can be matched with GLY and VAL but does not imply that GLY 
and VAL can be matched together. Using this syntax you can specify, for each 
residue, the allowed substitutions instead of simply defining groups of residues
which are all equivalent.

When you specify, as in the above examples, an equivalence that involves the
whole residue (i.e. you have omitted the .[atom] part), the points that make
up each residue will be paired in the same order as they appear in the
definition. e.g.

	def ALA=CA;O
	def VAL=CA;O

	equiv ALA=VAL

means that, in superimposing ALA and VAL, the program will pair CA with CA and
O with O. If you define the residue as follows

	def ALA=CA;O
	def VAL=O;CA

	equiv ALA=VAL

The CA and O of ALA will be paired with the O and CA of VAL, respectively.

If you want you can specify equivalences between specific residue fragments as
defined in the residue definition section. For instance if your residue definition
looks like this

	def=ALA=CA;avg(side_chain):bar
	def ASP=CA;CB;OG1

you could write the following equivalence

	equiv ALA.bar=ASP.OG1

You can also specify multiple atoms together, e.g.:

	equiv ALA.CA-bar=ASP.CA-OG1

Once again atoms will be paired according to the order in which they are written
(i.e. in the above case CA with CA and bar with OG1).
The atom names that you can use in these statements are the ones that you have
defined when specifying protein residues (see above).

- Wildcards

In specifying residue definitions and equivalences you can use wildcards. Two
types of wildcards are available "*" and "\". "*" is used in residue names and
it means "Any residue that does not match a more specific definition". "\" is
used in atom names with the meaning "Any atom whose name contains this
string". We will now describe in more detail how this wildcards apply to residue
definitions and equivalences. Please note that the only wildcard that can be
used in "res" mode is "*" when specifying residue definitions. Conversely all
the wildcards can be used in "res_multiple" and "atom" modes (see below).

The "*" in residue definitions is used to provide a definition for all the
remaining residues for which a more specific definition is not present. For
instance 

	def ALA=CA;CB
	def ASP=CA;avg(side_chain):bar
	def *=CA

means that ALA and ASP will be represented as specified, while all the other
residues will be defined by the CA only. Note that it is always a good idea to
include a catch all definition (i.e. using "*") at the end of the more
specific statements because the program will exit when a residue, for which no
definition is available, is encountered while parsing the PDB file.

The "\" in atom names is expanded, by looking at the PDB residue dictionary
(see "Command line syntax and files needed by superpose"), to all the atoms
whose name matches the string following the "\". For instance

	def ASP=\O

is equivalent to
	
	def ASP=O;OD1;OD2

given that the line in the PDB residue dictionary for ASP is

	ASP     ASP     D       N,CA,C,O,CB,CG,OD1,OD2

The "\" inside an "avg" clause are treated in the same way.

Similarly, in specifying equivalences, the "*" means "Any residue for which a
specific equivalence has not been provided" and "\" "All the atoms whose name
matches the string following the "\". 

These wildcards are especially useful if one wants to do atom-based
comparisons. For instance the following parameters file

	mode atom

	def *=\N;\O

	equiv *.\N=*.\N
	equiv *.\O=*.\O

means that all the residues should be represented with their nitrogen and
oxygen atoms and that only atoms of the same type are allowed to match.

- Modes

superpose can operate in three different modes which are explained below:

- res: In this mode the program will always output matches between pairs of
  residues. You'll never see cases where the same residue matches two
  different residues with different pseudoatoms. Moreover In the equivalences 
  section you can only specify one way to superimpose the same pair of residues.
  In other words something like this:

  	equiv ALA.CA=ASP.CA
	equiv ALA.CB=ASP.CB
	
  would be invalid because you have specified two alternative ways to
  superimpose the same pair of residues, ALA and ASP.

- res_multiple: this is similar to res mode. The difference is that in this
  mode you can specify multiple ways to superimpose the same pair of residues.
  During the computation the algorithm will select, for each pair of
  residues of the structure, what is the optimal way to superimpose them.

- atom: this is similar to res_multiple. The difference is that in this mode
  all the pseudoatoms that define a residue are treated independently.
  Therefore the program may output matches where the same residue is paired
  with two or more different residues using different pseudoatoms. For
  instance you might have the CA of res 1 matching with the CA of res 2 in the
  other structure and the CB of res 1 matching with the CB of res 3. In other
  words in this mode the notion of residue is abolished and every pseudoatom
  you have defined is treated as a different entity.

Aside from these differences the additional thing to consider is that mode
"res" is slightly faster than the other two. Modes are specified in the
parameter file (See "Specifying options and parameters").

- Output

The output is structured as follows:

	[match_id] [probe] [taget] [probe size] [target size] [rmsd] [length] [residue pairs] [trans]

The match_id is just a numeric id for the match, each match is located on a
different line. Probe and target are the chains involved in the match. Probe
and target sizes are the number of elements that can be independently matched
comprising the probe and the target chains. This is equal to the number of
residues for the "res" and "res_multiple" modes and the total number of
pseudoatoms in "atom" mode. This is useful if you want to calculate statistics
based e.g. on the ratio between the size of the match and the size of the
chains used as input.
The residue pairs give the correspondences between the two structures.  If you
are in "res_multiple" or "atom" and you have specified equivalences between
fragments of residues the name of the aminoacid will be followed by the name
of the fragment that has matched.  Trans is a set of 15 numbers which
represent the geometric transformation to be applied to the original pdb files
in order to superimpose the structures according to the match. The first 3
numbers are the translation along the x, y and z axis to be applied to the
probe. Similarly the following 3 represent the translation for the target. The
last 9 numbers are the three rows of the 3x3 rotation matrix which has to be
applied, after the translation, to the target _only_.

- Superimposing PDB files and visualising the structural similarities

Superpose will compare the two protein structures but it will not create
coordinate files with the structures roto-translated according to the
structural match. In order to transform the PDB files you need the program
superpose_pdb_files which is built together with superpose and should be
located in the build/src subdir (see, "Installation"). The syntax of this
program is as follows:

	superpose_pdb_files [-i input_file] [-d list of pdb file dir]

input_file is the output file from superpose, containing the matches that you
want to visualize, "list of pdb file dir" is a list of directories containing
the original pdb files which have been used as input to superpose. The program
will look in each directory trying to find the files and abort with an error
if it fails.  For each line in the input file (i.e. each match) this program
will create three files:

	[probe]_rot[match_id].ent
	[target]_rot[match_id].ent
	chim_[match_id].cmd

The first two files are the PDB files of the probe and target rotated and
translated according to the structural match having the given match_id. You
can visualize them using any molecular graphics program of your choice.
The third file is a script to be used with the excellent molecular viewer
Chimera (available for free at http://www.cgl.ucsf.edu/chimera/). To run the
script in chimera place all the three files in the same directory, open
chimera and fire up the command line (menu favorites -> command line)
In the command line type

	cd [name of dir where you put the files]
	source [name of the file containing the chimera script]

You should see the two structures with the residues comprising the match
represented as sticks and the remainder of the proteins depicted in ribbon
representation.







