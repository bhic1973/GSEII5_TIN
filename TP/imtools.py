import os
import glob
import random
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter


class ImageDegrader:
    def __init__(self, root_folder, rank, extensions=["*.jpg", "*.png", "*.jpeg"]):
        """
        Initialize the ImageDegrader class.

        Args:
            root_folder (str): Root folder to search for images.
            extensions (list): List of image file extensions to search for.
        """
        self.root_folder = root_folder
        self.file_exts = extensions
        self.rank = rank
        self.image_paths = self._collect_image_paths()
        self.image = self.choose_random_image()

    def _collect_image_paths(self):
        """
        Collect image file paths from the root folder and subfolders.

        Args:
            root_folder (str): Folder to search for images.
            extensions (list): List of file extensions to include.

        Returns:
            list: List of image file paths.
        """
        image_paths = []
        for ext in self.file_exts:
            image_paths.extend(glob.glob(os.path.join(self.root_folder, "**", ext), recursive=True))
        return image_paths

    def choose_random_image(self):
        """
        Choose a random image from the collected paths.

        Returns:
            PIL.Image.Image: A randomly chosen image.
        """
        rank = self.rank
        random.seed(rank+321)
        if not self.image_paths:
            raise ValueError("No image paths found.")
        random_path = random.choice(self.image_paths)
        return Image.open(random_path)

    def add_gaussian_noise(self, image, mean=-1):
        """
        Add Gaussian noise to an image.

        Args:
            image (PIL.Image.Image): Input image.
            mean (float): Mean of the Gaussian noise.
            std (float): Standard deviation of the Gaussian noise.

        Returns:
            PIL.Image.Image: Image with added noise.
        """
        random.seed()
        std = random.randint(10,50)
        im_hsv = image.convert('HSV')
        H, S, V = im_hsv.split()
        img_array = np.array(V,dtype=np.float32)
        noise = np.random.normal(mean, std, img_array.shape)
        noisy_array = np.clip(img_array + noise, 0, 255).astype(np.uint8)
        return Image.merge('HSV',[H,S,Image.fromarray(noisy_array)]).convert('RGB')

    def add_impulse_noise(self, image, noise_type="s&p", amount=0.05):
        """
        Add impulse noise (s, p, or s&p) to an image.

        Args:
            image (PIL.Image.Image): Input image.
            noise_type (str): Type of noise ("salt", "pepper", or "salt_and_pepper").
            amount (float): Proportion of pixels to be affected by noise (0 to 1).

        Returns:
            PIL.Image.Image: Image with added impulse noise.
        """
        # Convert the image to a NumPy array
        im_hsv = image.convert('HSV')
        H, S, V = im_hsv.split()
        img_array = np.array(V)
        random.seed()
        # Determine the number of pixels to modify
        total_pixels = img_array.size
        num_noisy_pixels = int(amount * total_pixels)

        # Generate random pixel indices
        row_indices = np.random.randint(0, img_array.shape[0], num_noisy_pixels)
        col_indices = np.random.randint(0, img_array.shape[1], num_noisy_pixels)

        # Add noise
        if noise_type == "s":  # Add only white pixels
            img_array[row_indices, col_indices] = 255
        elif noise_type == "p":  # Add only black pixels
            img_array[row_indices, col_indices] = 0
        elif noise_type == "s&p":  # Add both white and black pixels
            for i in range(num_noisy_pixels):
                if np.random.random() < 0.5:
                    img_array[row_indices[i], col_indices[i]] = 255  # Salt
                else:
                    img_array[row_indices[i], col_indices[i]] = 0  # Pepper

        # Convert back to a PIL image
        im_hsv = Image.merge('HSV',[H,S,Image.fromarray(img_array)])
        return im_hsv.convert('RGB')
 
    
    def apply_blur(self, image, sigma=1):
        """
        Apply Gaussian blur to an image.

        Args:
            image (PIL.Image.Image): Input image.
            sigma (float): Standard deviation for Gaussian kernel.

        Returns:
            PIL.Image.Image: Blurred image.
        """
        # Convert the image to a NumPy array
        im_hsv = image.convert('HSV')
        H, S, V = im_hsv.split()
        img_array = np.array(V,dtype=np.float32)
        # img_array = np.array(image, dtype=np.float32)
        blurred = gaussian_filter(img_array, sigma=sigma)
        blurred = np.clip(blurred,0,255).astype('uint8')
        im_hsv = Image.merge('HSV',[H,S,Image.fromarray(blurred)])
        return im_hsv.convert('RGB')

    def degrade_illumination(self, image, gamma=0.5, shift=-50):
        """
        Degrade illumination by applying gamma transformation and histogram shifting.

        Args:
            image (PIL.Image.Image): Input image.
            gamma (float): Gamma correction factor (>0 darkens the image, <1 lightens it).
            shift (int): Value to shift pixel intensity (-ve to darken, +ve to lighten).

        Returns:
            PIL.Image.Image: Image with degraded illumination.
        """
        # Convert the image to a NumPy array
        im_hsv = image.convert('HSV')
        H, S, V = im_hsv.split()
        img_array = np.array(V,dtype=np.float32)

        # Normalize the image to [-1, 1] for gamma correction
        img_array = img_array / 255.0
        img_array = np.clip(img_array, 0, 1)  # Ensure no overflow after normalization

        # Apply gamma correction
        gamma_corrected = np.power(img_array, gamma)

        # Scale back to [-1, 255] and apply histogram shift
        shifted = gamma_corrected * 255 + shift
        shifted = np.clip(shifted, 0, 255).astype(np.uint8)  # Ensure valid range

        # Convert back to a PIL image
        im_hsv = Image.merge('HSV',[H,S,Image.fromarray(shifted)])
        return im_hsv.convert('RGB')

    def apply_all_degradations(self):
        """
        Apply all types of degradation (noise, blur, illumination degradation) to the image.

        Args:
            image (PIL.Image.Image): Input image.

        Returns:
            PIL.Image.Image: Degraded image with all degradations applied sequentially.
        """

        # Apply blur
        degraded_image = self.apply_blur(self.image, sigma=random.uniform(0, 3))

        # Add noise
        noise_type = random.choice(['gaussian', 'impulsional'])
        if noise_type == 'gaussian':
            degraded_image = self.add_gaussian_noise(degraded_image)
        else:
            degraded_image = self.add_impulse_noise(degraded_image,noise_type=random.choice(['s&p', 's','p']), amount=random.uniform(0.05,0.1))

        # Degrade illumination
        #random.seed()
        degraded_image = self.degrade_illumination(degraded_image)

        return degraded_image
    
