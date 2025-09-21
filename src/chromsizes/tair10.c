/*
 *   bwtk: bigWig Toolkit
 *   Copyright (C) 2025  Benjamin Jean-Marie Tremblay
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

// To add additional preset genomes, simply copy this file into a new one
// inside src/chromsizes/, modify the values, and modify the chromSizesFromPreset()
// function inside src/bwtk.c to check for user input & #include the new file.

chromSizes->n = 7;
initChromSizes(chromSizes, chromSizes->n);
if (!ucsc) {
    const char *chroms[] = { "1", "2", "3", "4", "5", "Mt", "Pt" };
    copyStrings(chroms, chromSizes->chroms, chromSizes->n);
} else {
    const char *chroms[] = { "chr1", "chr2", "chr3", "chr4", "chr5", "chrM", "chrC" };
    copyStrings(chroms, chromSizes->chroms, chromSizes->n);
}
const uint32_t sizes[] = { 30427671, 19698289, 23459830, 18585056, 26975502, 366924, 154478 };
copyInts(sizes, chromSizes->sizes, chromSizes->n);

