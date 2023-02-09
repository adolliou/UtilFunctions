import os
from glob import glob


class MakeMovie:
    @staticmethod
    def encode(frame_folder: str, path_save: str, start_frame=0, nb_frame=None, suffix="*.png"):
        files = glob.glob(os.path.join(frame_folder, suffix))
        if nb_frame is None:
            nb_frame = len(files) - start_frame
        spawnline = 'ffmpeg -r 10 -s 640:480 -start_number "' + start_frame + '" -i "' + \
                    os.path.join(frame_folder, '%05d.png') \
                    + '" -vframes "' + nb_frame + '" -c:v libx264 -vf fps=25 -pix_fmt yuv420p -y "' + \
                    path_save + '"'
        os.system(spawnline)

