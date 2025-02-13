# Libraries imported
import matplotlib
matplotlib.use('tkAgg')
from matplotlib import pyplot as plt
from PIL import Image
import torch
from transformers import GLPNImageProcessor, GLPNForDepthEstimation

# 1. Getting feature extractore (GLPNImage) + model
# Code given:
feature_extractor = GLPNImageProcessor.from_pretrained("")
model = GLPNForDepthEstimation.from_pretrained("")