#!/bin/bash
# removesharedmemes.sh
# Author: Chris Carson
# Description: Removes .meme files from folder2 if same protein exists in folder1.
# Includes user confirmation before deletion.

if [ $# -ne 2 ]; then
    echo "Usage: $0 <folder1_path> <folder2_path>"
    echo "This will remove .meme files from folder2 if the same protein exists in folder1"
    echo "Example: $0 ./memes/set1 ./memes/set2"
    exit 1
fi

FOLDER1="$1"
FOLDER2="$2"

if [ ! -d "$FOLDER1" ]; then
    echo "Error: Folder '$FOLDER1' does not exist"
    exit 1
fi

if [ ! -d "$FOLDER2" ]; then
    echo "Error: Folder '$FOLDER2' does not exist"
    exit 1
fi

echo "Analyzing proteins in both folders..."

ls "$FOLDER1"/*.meme 2>/dev/null | xargs -n1 basename | cut -d'_' -f1 | sort -u > folder1_proteins.txt
ls "$FOLDER2"/*.meme 2>/dev/null | xargs -n1 basename | cut -d'_' -f1 | sort -u > folder2_proteins.txt

if [ ! -s folder1_proteins.txt ]; then
    echo "Warning: No .meme files found in $FOLDER1"
    rm -f folder1_proteins.txt folder2_proteins.txt
    exit 1
fi

if [ ! -s folder2_proteins.txt ]; then
    echo "Warning: No .meme files found in $FOLDER2"
    rm -f folder1_proteins.txt folder2_proteins.txt
    exit 1
fi

comm -12 folder1_proteins.txt folder2_proteins.txt > shared_proteins.txt

if [ ! -s shared_proteins.txt ]; then
    echo "No shared proteins found. Nothing to remove."
    rm -f folder1_proteins.txt folder2_proteins.txt shared_proteins.txt
    exit 0
fi

echo "Found $(wc -l < shared_proteins.txt) shared proteins:"
cat shared_proteins.txt

echo ""
echo "Files to be removed from $FOLDER2:"
echo "=================================="

while read -r protein; do
    files_found=$(ls "$FOLDER2"/${protein}_*.meme 2>/dev/null)
    if [ -n "$files_found" ]; then
        echo "$files_found" | while read -r file; do
            echo "$(basename "$file")"
        done
    fi
done < shared_proteins.txt

echo ""
read -p "Do you want to proceed with removing these files? (y/N): " confirm

if [[ $confirm =~ ^[Yy]$ ]]; then
    while read -r protein; do
        files_to_delete=$(ls "$FOLDER2"/${protein}_*.meme 2>/dev/null)
        if [ -n "$files_to_delete" ]; then
            echo "$files_to_delete" | while read -r file; do
                echo "Removing: $(basename "$file")"
                rm "$file"
            done
        fi
    done < shared_proteins.txt
    
    echo ""
    echo "Removal complete!"
    echo "Files removed from $(basename "$FOLDER2")"
else
    echo "Operation cancelled."
fi

rm -f folder1_proteins.txt folder2_proteins.txt shared_proteins.txt
