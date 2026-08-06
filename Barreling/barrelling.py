'''
Barrelling correction - G.Palombella - 06.24.2026
'''

import numpy as np
import cv2
import matplotlib.pyplot as plt
from scipy.cluster.vq import kmeans2
from optparse import OptionParser
import sys, os
import discorpy.proc.processing as proc
import discorpy.post.postprocessing as post
import discorpy.losa.loadersaver as losa

class barrelling:
    def __init__(self, options, grid_image): 
        self.options = options
        self.xcenter, self.ycenter, self.list_fact = self.compute_correction(grid_image)

    def red_mask(self, img, r_min=150, g_max=100, b_max=100):
        b, g, r = img[:, :, 0], img[:, :, 1], img[:, :, 2]
        mask = (r > r_min) & (g < g_max) & (b < b_max)
        return mask.astype(np.uint8)

    def assign_clusters_1d(self, vals, k, seed=0):
        # raggruppa coordinate 1D in k cluster -> assegna ogni punto alla sua riga/colonna
        vals = np.asarray(vals, dtype=float)
        init = np.percentile(vals, np.linspace(5, 95, k))  # k centri iniziali ben distribuiti
        centers, labels = kmeans2(vals, init, minit='matrix', seed=seed)
        order = np.argsort(centers)
        remap = {old: new for new, old in enumerate(order)}  # rinomina i cluster in ordine crescente
        labels_sorted = np.array([remap[l] for l in labels])
        return labels_sorted

    def extract_points(self, img, tol=15.0): 
        n_h, n_v = self.options.n_h, self.options.n_v
        mask = self.red_mask(img)
        _, _, _, centroids = cv2.connectedComponentsWithStats(mask, connectivity=8) # trova i centroidi delle componenti connesse --> i punti non sono puntiformi...
        pts = centroids[1:]

        xs, ys = pts[:, 0], pts[:, 1] 

        # conta vicini con y simile (riga orizzontale) e x simile (colonna verticale)
        count_h = np.array([np.sum(np.abs(ys - y) < tol) for y in ys]) # prende un punto di ys e misura la sua distanza da tutti gli altri punti, se questa è minore di una soglia, segno True. poi sommo tutti i True (1) e ottengo il valore per quel punto di ys, poi ripeto.
        count_v = np.array([np.sum(np.abs(xs - x) < tol) for x in xs])
        is_horizontal = count_h > count_v
        h_pts = pts[is_horizontal]
        v_pts = pts[~is_horizontal]

        h_line_labels = self.assign_clusters_1d(h_pts[:, 1], n_h)
        v_line_labels = self.assign_clusters_1d(v_pts[:, 0], n_v)

        h_xs = [h_pts[h_line_labels == i, 0].tolist() for i in range(n_h)]
        h_ys = [h_pts[h_line_labels == i, 1].tolist() for i in range(n_h)]
        v_xs = [v_pts[v_line_labels == i, 0].tolist() for i in range(n_v)]
        v_ys = [v_pts[v_line_labels == i, 1].tolist() for i in range(n_v)]

        return h_xs, h_ys, v_xs, v_ys

    def to_discorpy_format(self, xs_list, ys_list, sort_axis=0): 
        out = []
        for xs, ys in zip(xs_list, ys_list):
            arr = np.stack([np.asarray(ys, dtype=float), np.asarray(xs, dtype=float)], axis=1)
            if arr.shape[0] > 0:
                arr = arr[np.argsort(arr[:, sort_axis])]  # ordina i punti lungo la linea
                out.append(arr)
        return out

    def compute_center(self, img, list_hor_lines, list_ver_lines):
        H, W = img.shape[:2]
        xcenter0, ycenter0 = W / 2.0, H / 2.0  # stima iniziale: centro geometrico
        dists = []
        for line in list_hor_lines:
            if line.shape[0] > 1:
                dists.extend(np.abs(np.diff(line[:, 1])))  # diff sulle x (punti ordinati per x)
        point_dist = float(np.median(dists)) if dists else 50.0

        xcenter, ycenter = proc.find_cod_fine(list_hor_lines, list_ver_lines, xcenter0, ycenter0, point_dist)
        
        return xcenter, ycenter 

    def compute_correction(self, img, verbose=True): 
        h_xs, h_ys, v_xs, v_ys = self.extract_points(img)
        list_hor_lines = self.to_discorpy_format(h_xs, h_ys, sort_axis=1)
        list_ver_lines = self.to_discorpy_format(v_xs, v_ys, sort_axis=0)  
        xcenter, ycenter = self.compute_center(img, list_hor_lines, list_ver_lines)

        list_fact = proc.calc_coef_backward(list_hor_lines, list_ver_lines, xcenter, ycenter, 5) # 5 è il numero di coefficienti

        if verbose:
            print(f"Center of distortion: ({xcenter:.3f}, {ycenter:.3f})")
            print(f"Coefficients: {list_fact}")

        return xcenter, ycenter, list_fact

if __name__ == '__main__':
    parser = OptionParser()
    (options, args) = parser.parse_args()

    f = open(args[0], "r")
    params = eval(f.read())
    #f.close()

    for k, v in params.items():
        setattr(options, k, v)

    grid_image = cv2.imread(options.grid_image)  
    H_grid, W_grid = grid_image.shape[:2]

    os.makedirs(options.output_dir, exist_ok=True)

    out = barrelling(options, grid_image)

    losa.save_metadata_txt(f'{options.output_dir}/correction_coeffs.txt', out.xcenter, out.ycenter, out.list_fact)

    if options.correct:
        def unwarp(image):
            return post.unwarp_image_backward(image.astype(float), xcenter_t, ycenter_t, out.list_fact_t)
        
        to_correct = np.load(options.to_correct)
        H_target, W_target = to_correct.shape[:2]
        scale_x = W_target / W_grid
        scale_y = H_target / H_grid
        xcenter_t = out.xcenter * scale_x
        ycenter_t = out.ycenter * scale_y
        list_fact_t = [f / scale_x**i for i, f in enumerate(out.list_fact)]
        
        if getattr(options, 'ped', None):
            ped = np.load(options.ped)
            to_correct = to_correct - ped  # subtract pedestal
        corrected = unwarp(to_correct)
        
        np.save(f'{options.output_dir}/corrected.npy', corrected)
        plt.imsave(f'{options.output_dir}/corrected.png', corrected, cmap='gray')

    sys.exit(0)
