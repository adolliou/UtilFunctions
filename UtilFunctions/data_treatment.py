import copy
import cv2
import numpy as np


class BackgroundEstimation:
    @staticmethod
    def inpaint(data: np.array, mask: np.array):
        bound = copy.deepcopy(data)

        mini = bound.min()
        bound -= mini
        bound **= 0.5
        maxi = bound.max()
        bound /= maxi
        bound *= 255
        bound = bound.astype(np.uint8)
        mask = mask.astype(np.uint8)
        cv2.inpaint(bound, mask, 3, cv2.INPAINT_NS, bound)  # input and ouput image : bound
        bound = bound.astype(data.dtype)
        bound /= 255
        bound *= maxi
        bound **= 2
        bound += mini
        inpainted = np.copy(bound)  # remise aux bonnes valeurs de bound
        inpainted[inpainted < 0] = 0
        return inpainted
