#!/usr/bin/env python3
import os
import sys
import csv
import math
import glob
import numpy as np
from typing import List, Tuple

#import tensorflow as tf
from tensorflow.keras.models import load_model
#from tensorflow.keras.utils import load_img, img_to_array
from tensorflow.keras.preprocessing.image import load_img, img_to_array

# -----------------------------
# Config — edit as needed
# -----------------------------
MODEL_PATH = "image_classifier_model.h5"   # path to your saved model
IMAGES_DIR = "/path/to/new_images"         # folder containing images (flat or nested)
OUTPUT_CSV = "predictions.csv"

IMG_HEIGHT, IMG_WIDTH = 128, 128           # must match training
BATCH_SIZE = 32
THRESHOLD = 0.005                            # 0.005 => positive (class 1)
VALID_EXTS = {".jpg", ".jpeg", ".png", ".bmp", ".tif", ".tiff"}


# -----------------------------
# (Optional) make TF GPU friendlier
# -----------------------------
# def _setup_gpu_memory_growth():
#     try:
#         gpus = tf.config.list_physical_devices("GPU")
#         for gpu in gpus:
#             tf.config.experimental.set_memory_growth(gpu, True)
#     except Exception:
#         pass


# -----------------------------
# Utils
# -----------------------------
def list_image_paths(root: str, recursive: bool = True) -> List[str]:
    patterns = ["**/*"] if recursive else ["*"]
    all_paths = []
    for pat in patterns:
        all_paths.extend(glob.glob(os.path.join(root, pat), recursive=recursive))
    img_paths = [p for p in all_paths if os.path.splitext(p)[1].lower() in VALID_EXTS and os.path.isfile(p)]
    img_paths.sort()
    return img_paths


def preprocess_image(path: str, target_size: Tuple[int, int]) -> np.ndarray:
    # load_img returns PIL image; ensure RGB
    img = load_img(path, target_size=target_size)  # default color_mode='rgb'
    arr = img_to_array(img)                         # float32 array in [0..255]
    arr = arr / 255.0                               # rescale like training
    return arr


def batch(iterable, n=1):
    length = len(iterable)
    for ndx in range(0, length, n):
        yield iterable[ndx : ndx + n]


# -----------------------------
# Prediction
# -----------------------------
def main(model_path: str, images_dir: str, out_csv: str):
    if not os.path.exists(model_path):
        print(f"Model not found: {model_path}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isdir(images_dir):
        print(f"Images directory not found: {images_dir}", file=sys.stderr)
        sys.exit(1)

    #_setup_gpu_memory_growth()
    print("Loading model…")
    model = load_model(model_path)

    print("Scanning images…")
    paths = list_image_paths(images_dir, recursive=True)
    if not paths:
        print(f"No images found under {images_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(paths)} images. Predicting in batches of {BATCH_SIZE}…")
    probs = np.zeros((len(paths),), dtype=np.float32)

    idx = 0
    for chunk in batch(paths, BATCH_SIZE):
        # Preprocess this batch
        batch_arr = np.stack([preprocess_image(p, (IMG_HEIGHT, IMG_WIDTH)) for p in chunk], axis=0)
        # Predict: model outputs sigmoid for class 1 (positive)
        batch_probs = model.predict(batch_arr, verbose=0).reshape(-1)
        probs[idx : idx + len(chunk)] = batch_probs
        idx += len(chunk)

    # Convert to labels using threshold
    labels = np.where(probs >= THRESHOLD, "positive", "negative")

    # Write CSV
    print(f"Writing {out_csv}…")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["filepath", "prob_positive", "predicted_label"])
        for p, pr, lb in zip(paths, probs, labels):
            w.writerow([p, f"{pr:.6f}", lb])

    # Console preview (first 10)
    print("\nSample predictions:")
    for p, pr, lb in list(zip(paths, probs, labels))[:10]:
        print(f"{lb:8s}  prob={pr:.3f}  {p}")

    # Quick summary
    n_pos = int((labels == "positive").sum())
    n_neg = len(labels) - n_pos
    print(f"\nSummary: {n_pos} positive, {n_neg} negative (threshold={THRESHOLD})")
    print(f"Done. Results saved to: {out_csv}")


if __name__ == "__main__":
    # Allow overriding via CLI: python predict_images.py /path/to/model.h5 /path/to/images out.csv
    if len(sys.argv) >= 2:
        MODEL_PATH = sys.argv[1]
    if len(sys.argv) >= 3:
        IMAGES_DIR = sys.argv[2]
    if len(sys.argv) >= 4:
        OUTPUT_CSV = sys.argv[3]

    main(MODEL_PATH, IMAGES_DIR, OUTPUT_CSV)

#conda activate py38
#python predict_img.py ./denovo_MEI/train_model/all_merged/image_classifier_model.h5 ./denovo_MEI/evaluate_model/TPs prediction_tp.csv
