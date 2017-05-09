
# Usage
Source files are in `src/`

## Modified images
A database image can be modified in three ways : `imgRotate`, `imgScale`, `imgAddNoise`
To test the classification of randomly modified images, type `zsh scriptClassification.sh`. The recognition works when Rank is 1. Otherwise, we should do better.

## Distances
To show the `9` neighbors of a given image like `beetle-1`, do `python3 src/distance.py database/beetle-1.pgm 9 v`
Removing `v` gives the probability of matching for each class of `classes.csv`


# Similarity measure between two images

- values must be in [0,1]
- 1 means "similar"


# Classification evaluation

- the tool must output (standard output) a list of 70 values between 0
and 1 (classification score or probability to belong to the class)

- ```getRank``` returns the rank of a given class in a given
  classification score.
