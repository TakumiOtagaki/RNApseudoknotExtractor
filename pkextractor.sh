# module load rnaview
# export PDB_ID="1KPD"
# to call this script:
    # ./pkextractor.sh --pdb_id 1KPD(=PDB_ID) --chimera_or_pyMOL(=viewer) chimera
export PDB_ID=$1

# activate python venv
source /home/otgk/PKExtractor/venv/bin/activate



# argparser
if [ $# -ne 2 ]; then
    echo "Usage: $0 PDB_ID viewer"
    echo "viewer must be chimera or pyMOL"
    exit 1
fi

if [ $2 = "chimera" ]; then
    viewer="chimera"
elif [ $2 = "pyMOL" ]; then
    viewer="pyMOL"
else
    echo "viewer must be chimera or pyMOL"
    exit 1
fi

if [ $viewer == "chimera" ]; then
    export def="colordef"
elif [ $viewer == "pyMOL" ]; then
    export def="set_color"
fi



echo "PDB_ID:$PDB_ID"

# at pdb directory
cd /home/otgk/PKExtractor/PKExtractor/pdb
echo "fetching ${PDB_ID}.pdb"
rna_pdb_tools.py --fetch $PDB_ID

echo "${PDB_ID}.pdb downloaded"

cd ..
echo "get sequence from ${PDB_ID}.pdb"
rna_pdb_tools.py --get-seq pdb/$PDB_ID.pdb > seq/${PDB_ID}.fasta
echo "sequence extraction done"

# if seq/${PDB_ID}.fasta is empty: exit
if [ -z "$(cat seq/${PDB_ID}.fasta)" ]; then
    echo "seq/${PDB_ID}.fasta is empty. exit"
    echo "This may be because ${PDB_ID}.pdb is not RNA structure."
    exit 1
fi

echo "get canonical basepair information from ${PDB_ID}.pdb by using rnaview"
cd pdb
rnaview ${PDB_ID}.pdb
echo "rnaview process done"

cd ..
# cat ${PDB_ID}.pdb.out | grep A: | sed -e "s/ \+/\t/g" | awk -F "\t" '(($7="W/W") && ($8="cis")) || ($9='XIX') || ($9="XXVIII") {print $0}' | sed -e "s/ /\t/g" | cut -f 2,5 | sed -e "s/,//g" | sed -e "s/_/\t/g" | sed -e "s/-/\t/g" > ../CanonicalBP/${PDB_ID}.tsv
# cd .. # move to PKExtractor directory
cat pdb/${PDB_ID}.pdb.out | awk '/BEGIN_base-pair/,/END_base-pair/' |  head -n -1 | tail -n +2 | sed -e "s/ \+/\t/g"  | sed -e "s/^\t//g" | sed -e "s/://g" | 
 awk -F "\t" '(($2==$6)) {print $0}' | # chain-internal bp
 awk -F "\t" '(($9=="XX") || ($9=="XIX") || ($9=="XXVIII")) {print $0}' > CanonicalBP/${PDB_ID}.tsv # canonical bp


# for each chain, extract canonical bp and write to file. then, PseudoKnot extractor will read this file.
cat seq/${PDB_ID}.fasta | grep ">" | cut -d ":" -f 1 | sed -e "s/>//g" |
while read chain
do
    echo "chain:$chain is processing"
    cat CanonicalBP/${PDB_ID}.tsv | awk -F "\t" -v chain=$chain '(($2==chain) && ($6==chain)) {print $3,$5}' > CanonicalBP/${PDB_ID}_${chain}.tsv
    # PKExtractor.py ../CanonicalBP/${PDB_ID}_${chain}.tsv ../
    python PKExtractor.py --input_file CanonicalBP/${PDB_ID}_${chain}.tsv --output_file PKresult/${PDB_ID}_${chain}.tsv

    # generate chimera script. coloring by PK id (column 3)
    echo "open ../PKresult/${PDB_ID}_${chain}.tsv and generate chimera script: coloring/${viewer}_${PDB_ID}_${chain}.txt"
    # if chimera: open pdb elif pyMOL: fetch pdb
    if [ $viewer = "chimera" ]; then
        echo "open ${PDB_ID}" > coloring/${viewer}_${PDB_ID}_${chain}.txt
    elif [ $viewer = "pyMOL" ]; then
        echo "fetch ${PDB_ID}" > coloring/pyMOL_${PDB_ID}_${chain}.txt
    fi
    # echo "open ${PDB_ID}" > coloring/${viewer}_${PDB_ID}_${chain}.txt
    
    # if chimera: colordef elif pyMOL: set_color. (color definion)


    echo "# define colors for each layer" >> coloring/${viewer}_${PDB_ID}_${chain}.txt
    echo "$def NonBP white" >> coloring/${viewer}_${PDB_ID}_${chain}.txt
    echo "$def layer0 gray" >> coloring/${viewer}_${PDB_ID}_${chain}.txt
    echo "$def layer1 red" >> coloring/${viewer}_${PDB_ID}_${chain}.txt
    echo "$def layer2 blue" >> coloring/${viewer}_${PDB_ID}_${chain}.txt
    echo "$def layer3 orange" >> coloring/${viewer}_${PDB_ID}_${chain}.txt
    echo "$def layer4 green" >> coloring/${viewer}_${PDB_ID}_${chain}.txt
    echo "$def layer5 purple" >> coloring/${viewer}_${PDB_ID}_${chain}.txt
    echo "$def layer6 magenta" >> coloring/${viewer}_${PDB_ID}_${chain}.txt
    echo "$def layer7 yellow" >> coloring/${viewer}_${PDB_ID}_${chain}.txt
    echo "$def layer8 cyan" >> coloring/${viewer}_${PDB_ID}_${chain}.txt
    echo "$def layer9 brown" >> coloring/${viewer}_${PDB_ID}_${chain}.txt



    # coloring
    echo "# color by PK id" >> coloring/${viewer}_${PDB_ID}_${chain}.txt
    echo "color NonBP"  >> coloring/${viewer}_${PDB_ID}_${chain}.txt # whilte.
    for i in {0..9}
    do
        
        # if (cat PKresult/${PDB_ID}_${chain}.tsv | awk -F "\t" -v i=$i '(($3==i)) {print $1,$2}' is empty): skip
        if [ -z "$(cat PKresult/${PDB_ID}_${chain}.tsv | awk -F "\t" -v i=$i '(($3==i)) {print $1,$2}')" ]; then
            continue
        fi

        if [ $viewer == "chimera" ]; then
            # for command in ["color", "bondcolor"]
            for command in "color" "bondcolor"
            do
                echo -n "$command layer${i} " >> coloring/${viewer}_${PDB_ID}_${chain}.txt
                cat PKresult/${PDB_ID}_${chain}.tsv | awk -F "\t" -v i=$i '(($3==i)) {print $1,$2}' |  sed -e "s/ /\n/g" |
                while read line
                do
                    echo -n ":$line.$chain " >> coloring/${viewer}_${PDB_ID}_${chain}.txt
                done
                echo "" >> coloring/${viewer}_${PDB_ID}_${chain}.txt
            done
        elif [ $viewer == "pyMOL" ]; then
            echo -n "color layer${i}, " >> coloring/${viewer}_${PDB_ID}_${chain}.txt
            # color layer$i, (resi i and chain $chain)
            cat PKresult/${PDB_ID}_${chain}.tsv | awk -F "\t" -v i=$i '(($3==i)) {print $1,$2}' |  sed -e "s/ /\n/g" |
            while read line
            do
                echo "((resi $line) and chain $chain)" >> coloring/${viewer}_${PDB_ID}_${chain}.txt
            done
        fi
    done
    echo "generate chimera script done!"
    echo "show ./coloring/${viewer}_${PDB_ID}_${chain}.txt"
done


# concat all $viewer script
echo "concat all $viewer script"
if [ $viewer == "chimera" ]; then
    echo "open ${PDB_ID}" > coloring/${viewer}_${PDB_ID}.txt
    cat coloring/${viewer}_${PDB_ID}_A.txt | grep $def >> coloring/${viewer}_${PDB_ID}.txt
    cat coloring/${viewer}_${PDB_ID}_*.txt | grep -v "open ${PDB_ID}" | grep -v "color NonBP" | grep -v $def >> coloring/${viewer}_${PDB_ID}.txt
elif [ $viewer == "pyMOL" ]; then
    echo "fetch ${PDB_ID}" > coloring/${viewer}_${PDB_ID}.txt
    cat coloring/${viewer}_${PDB_ID}_A.txt | grep $def >> coloring/${viewer}_${PDB_ID}.txt
    cat coloring/${viewer}_${PDB_ID}_*.txt | grep -v "open ${PDB_ID}" | grep -v "color NonBP" | grep -v $def >> coloring/${viewer}_${PDB_ID}.txt
    cat coloring/${viewer}_${PDB_ID}_*.txt | grep -v "fetch ${PDB_ID}" | grep -v "color NonBP" >> coloring/${viewer}_${PDB_ID}.txt
fi

echo "all operation done!"
echo "show ./coloring/${viewer}_${PDB_ID}.txt"
deactivate
