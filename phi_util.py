from sunpy.visualization.colormaps.color_tables import hmi_mag_color_table, _cmap_from_rgb



class PHIUtil:
    @staticmethod
    
    def build_color_table_phi(pmax, pmin, threshold_small_scale, fs=0.1):
        """
        Multi-scale colorbar for HRT, adapted from hmimag colorbar. must be used with the right norm.

        0           A + fs*256           A = 127 - 3l*256     127     B = 127 + 3l*256        B + B + fs*256      pmax
        pmin        
        Args:
            pmax (_type_): _description_
            pmin (_type_): _description_
            sigma_b (_type_): _description_
        """    
        

        # path_csv = os.path.join(Path(__file__).parents[0], "hmi_mag.csv")
        # data_old  = np.loadtxt(path_csv)
        data_new = np.zeros((256, 3), dtype=int)

        l = 1/(pmax - pmin)
        A = int(127 - threshold_small_scale * l * 256)
        B = int(127 + threshold_small_scale * l * 256)
        factor = int(fs*256)
        assert factor < A
        points_new = [0,            A - factor-1,
                    A - factor,     A-1,
                    A,              B,
                    B+1,            B + factor - 1,
                    B + factor,     255]
        
        # data_new[]
        
        points_old = [0,            107,
                    108,            125,
                    126,            129,
                    130,            147,
                    148,            255]

        # values_old = [
        #     [137, 0, 20],           [255, 219, 188], 
        #     [255, 103, 0],          [255, 255, 0],  
        #     [199, 200, 130],        [95, 129, 104], 
        #     [0, 120, 14],           [0, 255, 92], 
        #     [200, 200, 255],        [0, 0, 73],   
        # ]

        values_old = [
                [0, 0, 73],             [200, 200, 255],
                [0, 255, 92],           [0, 120, 14],  
                [45, 45, 45],           [210, 210, 210],         
                [255, 255, 0],          [255, 103, 0], 
                [255, 219, 188],        [137, 0, 20],
            ]


        # breakpoint()
        for ll in range(5):
            n = points_new[2*ll+1] - points_new[2*ll] + 1
            # print(f"{n=}")
            for kk in range(3):
                data_new[points_new[2*ll]:(points_new[2*ll + 1] + 1), kk] = np.array(np.linspace(values_old[2*ll][kk], values_old[2*ll + 1][kk], n), dtype=int)
        return _cmap_from_rgb(data_new[:, 0], data_new[:, 1], data_new[:, 2], "hrtmap")




            