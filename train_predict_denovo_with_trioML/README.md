# xTea-trioML image-based post-filtering module for de novo retrotransposon calls

## Overview

This module adds an image-based machine learning post-filtering step on top of the **xTea-trioML** pipeline for de novo retrotransposon insertion discovery in trio data. The goal is to reduce false positives among candidate de novo calls by:

1. generating labeled training images from xTea-trioML trio outputs,
2. training a convolutional neural network (CNN) on those images, and
3. applying the trained model to new candidate screenshots to classify each candidate as **positive** or **negative**.

The current implementation supports the three TE classes used throughout the scripts:

- **L1**
- **SVA**
- **Alu**

In practice, the workflow starts from xTea-trioML result folders for proband, father, and mother. The `prep_denovo_training_data_from_trio.py` script prepares shell commands to generate screenshot images for both positive and negative training examples. Those images are then used by `train_model.py` to build a binary image classifier. Finally, `classify_images.py` applies the trained model to a folder of candidate images and outputs a CSV file with predicted labels and probabilities.

---

## Module contents

### `prep_denovo_training_data_from_trio.py`
Prepares training data generation commands from xTea-trioML trio outputs.

Main tasks performed by the script:

- loads trio-specific xTea output folders for **proband**, **paternal**, and **maternal** samples,
- loads a BAM list for all trios,
- loads a manually curated list of true-positive de novo insertions,
- selects high-quality heterozygous insertions in the proband as candidate **positive** examples,
- selects raw de novo calls from xTea output that are not in the manually curated TP set as candidate **negative** examples,
- creates `positive/` and `negative/` subdirectories under the working folder,
- generates shell scripts that can be executed later to produce screenshot images through the `XBamSnap` helper. 

Positive examples are chosen from high-confidence proband insertions after parent/proband comparisons and feature-based filtering. The code checks for overlap against parental calls and uses feature-derived clip/disc ratios to keep insertions whose VAF-like ratios fall within the configured range. It also attempts to pair each proband event with a pseudo-parental trio that does not contain an overlapping insertion, so the generated image mimics a de novo trio context. Negative examples are collected from `candidate_disc_filtered_cns_post_filtering_denovo.txt` after excluding manually inspected true positives. 

### `run_gnrt_screenshot_cmd.sh`
A demo shell script showing how to run `prep_denovo_training_data_from_trio.py`.

It defines example values for:

- the xTea code directory,
- working directory,
- reference genome,
- TE type,
- manually curated TP file,
- trio folder lists and BAM list.

It then runs the preparation script to emit a shell script named like:

- `generate_positive_training_images_L1.sh`

The preparation script also writes a companion negative-image script:

- `generate_positive_training_images_L1.sh.negative.sh`

### `train_model.py`
Trains a binary CNN image classifier using TensorFlow/Keras.

The current script:

- expects a directory structure with `train/positive/` and `train/negative/`,
- rescales images to `[0, 1]`,
- reserves 20% of the data for validation via `validation_split=0.2`,
- trains a small CNN with three convolution/max-pooling blocks,
- uses image size **128 x 128**,
- trains for **10 epochs**,
- saves the model to `image_classifier_model.h5`,
- reports validation accuracy and a classification report on the validation split.

### `classify_images.py`
Applies a trained model to a new folder of candidate images.

The script:

- loads a saved Keras model,
- recursively scans a directory for supported image files,
- resizes all images to **128 x 128**,
- rescales pixel values to `[0, 1]`,
- predicts a positive probability for each image,
- converts probabilities into binary labels using a user-editable threshold,
- writes a CSV file containing:
  - `filepath`
  - `prob_positive`
  - `predicted_label`

Supported input image extensions are:

- `.jpg`
- `.jpeg`
- `.png`
- `.bmp`
- `.tif`
- `.tiff`

---

## Expected inputs

### 1. Trio xTea output folder lists
The preparation script expects three text files describing xTea output folders:

- paternal folder list: `--pa`
- maternal folder list: `--ma`
- proband folder list: `--proband`

Each line should contain:

```text
<family_or_sample_id> <path_to_xtea_output_folder>
