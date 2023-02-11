import os
from glob import glob


class MakeMovie:
    @staticmethod
    def encode(frame_folder: str, path_save: str, start_frame=0, nb_frame=None, suffix="png", resolution="640:480",
               frames_ps=10):
        files = glob(os.path.join(frame_folder, f'*.{suffix}'))
        print(f'{files=}')
        if nb_frame is None:
            nb_frame = len(files) - start_frame
        # spawnline = f'ffmpeg -s {resolution} -start_number "' + f'{start_frame}' + '" -i "' + \
        #             os.path.join(frame_folder, f'%05d.{suffix}') \
        #             + '" -vframes "' + f'{nb_frame}' + '" -c:v libx264 "' + \
        #             path_save + '"'
        spawnline = f'ffmpeg -framerate 1 -s {resolution} -start_number "' + f'{start_frame}' + '" -i "' + \
                    os.path.join(frame_folder, f'%05d.{suffix}') \
                    + '" -vframes "' + f'{nb_frame}' + '" -c:v libx264 -r 1 -flvflags no_duration_filesize "' + \
                    path_save + '"'
        # spawnline = f'ffmpeg -r {frames_ps} -s {resolution} -start_number "' + f'{start_frame}' + '" -i "' + \
        #             os.path.join(frame_folder, f'%05d.{suffix}') \
        #             + '" -vframes "' + f'{nb_frame}' + '" -c:v libx264 -vf fps=25 -pix_fmt yuv420p -y "' + \
        #             path_save + '"'
        os.system(spawnline)

