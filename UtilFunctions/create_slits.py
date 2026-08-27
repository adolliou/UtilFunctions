import numpy as np
import astropy.units as u
import copy
from scipy.interpolate import CubicSpline



class CreateSlit:


    @staticmethod
    
    def extract_stack_plot_slit(unit, header_eui, p1_, p2_, width):
        if p1_[1] <= p2_[1]:
            p1 = p1_
            p2 = p2_
        else:
            p1 = p2_
            p2 = p1_

        x = np.array([p2[0] - p1[0], p2[1] - p1[1]])  # Vecteur p1 ==> p2
        if x[1] == 0:
            normal = np.array([0, 1])
        elif x[0] == 0:
            normal = np.array([1, 0])
        elif (x[0] == 0) & (x[1] == 0):
            raise NotImplementedError
        else:
            a = 1 / np.sqrt(1 + (x[0] / x[1]) ** 2)
            normal = np.array([a, - a * (x[0] / x[1])])  # vecteur normal a x

        dlong = u.Quantity(header_eui["CDELT1"], header_eui["CUNIT1"]).to(unit).value
        dlatg = u.Quantity(header_eui["CDELT2"], header_eui["CUNIT2"]).to(unit).value

        normx = np.sqrt(x[0] ** 2 + x[1] ** 2)
        dx = x / normx
        delt = np.mean([dlong, dlatg])
        ttx = None
        tty = None
        if p1[0] != p2[0]:
            ttx = np.arange(p1[0], p2[0], dx[0] * delt)
        if p1[1] != p2[1]:
            tty = np.arange(p1[1], p2[1], dx[1] * delt)
        if (ttx is not None) & (tty is not None):
            raise NotImplementedError

        if ttx is None:
            ttx = np.array([p1[0]] * len(tty))
        if tty is None:
            tty = np.array([p1[1]] * len(ttx))

        kk = np.arange(-0.5 * width, 0.5 * width, delt)
        ppx = np.zeros((len(kk), len(ttx)))
        ppy = np.zeros((len(kk), len(ttx)))

        for jj, el in enumerate(kk):
            ppx[jj, :] = ttx + el * normal[0]
            ppy[jj, :] = tty + el * normal[1]
        return ppx, ppy

    @staticmethod
    def extract_stack_plot_spline(
        unit                : str,
        header,
        width               : float,
        lon_points          : list,
        lat_points          : list,
        cut_lon             : list = None,
        cut_lat             : list = None, 
        angle               : float = None,
        shift_lonlat        : list = None,
        perp_float          : int = None
                                ):
        """Create slits from multiple successive splines. 

        Args:
            unit (_type_): unit output of the points
            header (_type_): fits header of the data where to create the slit. 
            width (_type_): with of the slit (in "unit")
            lon_points (list): list of coordinates points (in units) that will define the splines 
            lat_points (list): list of coordinates points (in units) that will define the splines 
            cut_lon (_type_, optional): size 2 lists being the boundaries the final slits [min, max] in longitude
            cut_lat (_type_, optional): size 2 lists being the boundaries the final slits [min, max] in latitude
            angle (_type_, optional): Either lon_points and lat_points should be strictly increasing. However, if this is not the case, 
                one can try rotating the slit with an angle so that one of the reprojected axis will be strictly increasing. 
                The slit is derotated after, so the angle parameter does not impact the slit coordinates.   
            shift_lonlat (_type_, optional): size 2 list [lon, lat] that shifts the slit by a given vector. 
            perp_float: choose an location along the slit [0.0 < x < 1.0]. Will return the slit perpendicular to the original
            slit at the index corresponding to this location. The width will be the same as the one provided as input.  

        Returns:
            ttx: (axis, vertical) 2D array of the slit points separated by the pixel size
        """
        do_rotate_angle = False
        lon_points = np.array(lon_points)
        lat_points = np.array(lat_points)

        if angle is not None:
            do_rotate_angle = True
            angle = np.deg2rad(angle)
            cos = np.cos(angle)
            sin = np.sin(angle)
            lon_points_ = copy.deepcopy(lon_points)
            lat_points_ = copy.deepcopy(lat_points)

            lon_points = lon_points_ * cos - lat_points_ * sin
            lat_points = lon_points_ * sin + lat_points_ * cos 


        is_increasing_lon = all(x < y for x, y in zip(lon_points, lon_points[1:]))
        is_increasing_lat = all(x < y for x, y in zip(lat_points, lat_points[1:]))

        if is_increasing_lon:
            cs = CubicSpline(lon_points, lat_points)
        elif is_increasing_lat:
            cs = CubicSpline(lat_points, lon_points)
        else:
            raise ValueError("longitude or latitude is not strictly increasing")


        dlong = u.Quantity(header["CDELT1"], header["CUNIT1"]).to(unit).value
        dlatg = u.Quantity(header["CDELT2"], header["CUNIT2"]).to(unit).value
        if is_increasing_lon:
            ttx_ = np.arange(lon_points[0], lon_points[-1] + 0.1 * dlong, 0.1 * dlong)
            tty_ = cs(ttx_)
        elif is_increasing_lat:
            tty_ = np.arange(lat_points[0], lat_points[-1] + 0.1 * dlatg, 0.1 * dlatg)
            ttx_ = cs(tty_)

        cdelt = np.mean([dlong, dlatg])
        dr = np.array([np.sqrt((ttx_[n + 1] - ttx_[n]) ** 2 + (tty_[n + 1] - tty_[n]) ** 2) for n in range(len(ttx_) - 1)])
        length_loop = dr.sum()
        dr_summed = [0]
        for dr_ in dr:
            dr_summed.append(dr_summed[-1] + dr_)
        # npoints_parallel = np.round(length_loop/cdelt) 
        dr_summed_regular = np.arange(0, length_loop + cdelt, cdelt)

        if is_increasing_lon:
            ttx = np.interp(x=dr_summed_regular, xp=dr_summed, fp=ttx_)
            tty = cs(ttx)

            dttx = ttx[-1] - ttx[-2]
            ttx_ext = np.append(ttx, ttx[-1] + dttx)
            tty_ext = cs(ttx_ext)
        elif is_increasing_lat:
            tty = np.interp(x=dr_summed_regular, xp=dr_summed, fp=tty_)
            ttx = cs(tty)

            dtty = tty[-1] - tty[-2]
            tty_ext = np.append(tty, tty[-1] + dtty)
            ttx_ext = cs(tty_ext)

        norm = np.zeros((len(ttx_ext) - 1, 2))
        for ii in range(len(ttx_ext) - 1):
            x = ttx_ext[ii + 1] - ttx_ext[ii]
            y = tty_ext[ii + 1] - tty_ext[ii]

            if is_increasing_lon:
                a = 1 / np.sqrt(1 + (x / y) ** 2)
                norm[ii, :] = np.array([a, - a * (x / y)])
            elif is_increasing_lat:
                a = 1 / np.sqrt(1 + (y / x) ** 2)
                norm[ii, :] = np.array([- a * (y / x), a ])

                
        ttx_plus = np.zeros_like(ttx)
        tty_plus = np.zeros_like(tty)
        ttx_minus = np.zeros_like(ttx)
        tty_minus = np.zeros_like(tty)

        ttx_tmp = np.array([ttx + 0.5 * width * norm[:, 0], ttx - 0.5 * width * norm[:, 0]])
        tty_tmp = np.array([tty + 0.5 * width * norm[:, 1], tty - 0.5 * width * norm[:, 1]])
        if is_increasing_lon:
            index_plus = tty_tmp.argmax(axis=0)
            index_minus = tty_tmp.argmin(axis=0)
        elif is_increasing_lat:
            index_plus = ttx_tmp.argmax(axis=0)
            index_minus = ttx_tmp.argmin(axis=0)

        for jj, idx_plus in enumerate(index_plus):
            ttx_plus[jj] = ttx_tmp[idx_plus, jj]
            tty_plus[jj] = tty_tmp[idx_plus, jj]
            idx_minus = index_minus[jj]

            ttx_minus[jj] = ttx_tmp[idx_minus, jj]
            tty_minus[jj] = tty_tmp[idx_minus, jj]


        len_ = np.sqrt((ttx_plus[0] - ttx_minus[0]) ** 2 + (tty_plus[0] - tty_minus[0]) ** 2)
        n_perp = int(len_ // cdelt)
        ttx = np.zeros((n_perp, len(ttx_plus)))
        tty = np.zeros((n_perp, len(tty_plus),))

        for jj in range(len(ttx_plus)):
            ttx[:, jj] = np.linspace(ttx_minus[jj], ttx_plus[jj], n_perp)
            tty[:, jj] = np.linspace(tty_minus[jj], tty_plus[jj], n_perp)



        if do_rotate_angle:
            do_rotate_angle = True
            angle_inv = -angle
            cos = np.cos(angle_inv)
            sin = np.sin(angle_inv)
            ttx_ = copy.deepcopy(ttx)
            tty_ = copy.deepcopy(tty)

            ttx = ttx_ * cos - tty_ * sin
            tty = ttx_ * sin + tty_ * cos         

        if cut_lon is not None:
            interval_lon = np.logical_and(np.nanmean(ttx, axis=0) >= cut_lon[0], np.nanmean(ttx, axis=0) <= cut_lon[1])
            ttx = ttx[:, interval_lon]
            tty = tty[:, interval_lon]        

        if cut_lat is not None:
            interval_lat = np.logical_and(np.nanmean(tty, axis=0) >= cut_lat[0], np.nanmean(tty, axis=0) <= cut_lat[1])
            ttx = ttx[:, interval_lat]
            tty = tty[:, interval_lat]  

        if shift_lonlat is not None:
            ttx += shift_lonlat[0]
            tty += shift_lonlat[1]

        if perp_float is not None:
            if (perp_float < 0) or (perp_float > 1):
                raise ValueError("perp_float input should be a float between 0.0 and 1.0")  
            len_along               = ttx.shape[1]
            index_perp              = np.round((len_along - 1) * perp_float)
            width_index             = ttx.shape[0]
            index_perp_lower        = index_perp - np.round(width_index/2)
            index_perp_upper        = index_perp + np.round(width_index/2)
            if (index_perp_lower < 0) or (index_perp_upper > len_along):
                raise ValueError("The perpendicular slit does not fit into the original slit : choose a lower width.")  

            # Cut the original slit so that the width of the final perpendicular is the same as the original slit.  
            ttx             = ttx[:, index_perp_lower:index_perp_upper]
            tty             = tty[:, index_perp_lower:index_perp_upper]

            # Transpose the matrix to create a perpendicular slit
            ttx             = ttx.T
            tty             = tty.T

        
        return ttx, tty
