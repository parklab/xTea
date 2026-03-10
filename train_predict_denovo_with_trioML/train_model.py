import os
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, accuracy_score
from tensorflow.keras.preprocessing.image import ImageDataGenerator
from tensorflow.keras import layers, models, optimizers

# Set these paths to your image folders
positive_dir = '/n/no_backup2/dbmi/park/simon_chu/denovo_MEI/train_model/all_merged/train/positive'
negative_dir = '/n/no_backup2/dbmi/park/simon_chu/denovo_MEI/train_model/all_merged/train/negative'
img_height, img_width = 128, 128  # Adjust as needed
batch_size = 32

# Prepare ImageDataGenerator for loading images
datagen = ImageDataGenerator(rescale=1./255, validation_split=0.2)

train_gen = datagen.flow_from_directory(
    directory=os.path.dirname(positive_dir),
    target_size=(img_height, img_width),
    batch_size=batch_size,
    classes=[os.path.basename(negative_dir), os.path.basename(positive_dir)],
    class_mode='binary',
    subset='training',
    shuffle=True
)

val_gen = datagen.flow_from_directory(
    directory=os.path.dirname(positive_dir),
    target_size=(img_height, img_width),
    batch_size=batch_size,
    classes=[os.path.basename(negative_dir), os.path.basename(positive_dir)],
    class_mode='binary',
    subset='validation',
    shuffle=False
)
print(os.path.basename(negative_dir), os.path.basename(positive_dir))
# Build a simple CNN model
model = models.Sequential([
    layers.Conv2D(32, (3,3), activation='relu', input_shape=(img_height, img_width, 3)),
    layers.MaxPooling2D(2,2),
    layers.Conv2D(64, (3,3), activation='relu'),
    layers.MaxPooling2D(2,2),
    layers.Conv2D(128, (3,3), activation='relu'),
    layers.MaxPooling2D(2,2),
    layers.Flatten(),
    layers.Dense(128, activation='relu'),
    layers.Dense(1, activation='sigmoid')
])

model.compile(optimizer=optimizers.Adam(learning_rate=1e-4),
              loss='binary_crossentropy',
              metrics=['accuracy'])

# Train the model
history = model.fit(
    train_gen,
    epochs=10,
    validation_data=val_gen
)

# Save the model
model.save('image_classifier_model.h5')

# Evaluate on validation set
val_gen.reset()
y_true = val_gen.classes
print(np.unique(y_true, return_counts=True))
y_pred = (model.predict(val_gen) > 0.5).astype(int).flatten()

print("Validation Accuracy:", accuracy_score(y_true, y_pred))
print(classification_report(y_true, y_pred, target_names=['negative', 'positive']))

