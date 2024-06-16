    # @property
    # def covariance_propagation_vector(self) -> np.ndarray:
    #     return self._covariance_propagation_vector

    # def propagate_covariance(self, covariance_propagation_vector: np.ndarray) -> None:
    #     self._covariance_propagation_vector = np.convolve(covariance_propagation_vector, self.covariance_propagation_vector)

    # def covariance(self, region: Region):
    #     propagation_hwidth = (self._covariance_propagation_vector.shape[0] - 1) // 2  # construct the propagation matrix using stride trick
    #     zero_pad = np.zeros(region.length - 1, dtype=self.covariance_propagation_vector.dtype)
    #     row_templet = np.concatenate((zero_pad, self.covariance_propagation_vector, zero_pad))
    #     stride = row_templet.strides[0]
    #     covariance_propagation_matrix = np.lib.stride_tricks.as_strided(row_templet[region.length - 1:], shape=(region.length, region.length + 2*propagation_hwidth), strides=(-stride, stride))
    #     left_width, right_width = min(region.left, propagation_hwidth), min(self.length - region.right - 1, propagation_hwidth)
    #     unpropagated_covariance = np.diag(self[region.left - left_width: region.right + right_width + 1])  # construct the unpropagated covariance matrix
    #     unpropagated_covariance = np.pad(unpropagated_covariance, (propagation_hwidth - left_width, propagation_hwidth - right_width), mode='edge')
    #     return np.einsum('ik,kg,jg->ij', covariance_propagation_matrix, unpropagated_covariance, covariance_propagation_matrix)

    # def slice(self, erg: tuple) -> "Spectrum":
    #     ergstart = erg[0] if erg[0] is not None else self.energies.min()
    #     ergstop = erg[1] if erg[1] is not None else self.energies.max()
    #     start = np.abs(self.energies - ergstart).argmin()
    #     stop = np.abs(self.energies - ergstop).argmin()
    #     return self[start:stop+1].view(Spectrum)

    # def copy(self):
    #     return deepcopy(self)

    # def estimate_area(self, region: Region) -> tuple:
    #     """
    #     Estimate peak area of spectrum in a region based on classic counting method and linear baseline assumption.
    #     """
    #     total = sum(self[region.indexes])
    #     len_left = min(2, region.left)
    #     len_right = min(2, self.length - region.right - 1)
    #     baseline = (self[region.left-len_left: region.left+1].mean() + self[region.right: region.right+len_right+1].mean()) * region.length / 2
    #     return total, baseline

    # def estimate_fitness(self, region: Region) -> float:
    #     """
    #     Estimate fitness of spectrum in a region based on fitted parameters.
    #     """
    #     fcounts = region.fit_peaks() + region.fit_baseline()
    #     return 1 - ((fcounts - self[region.indexes])**2).sum() / (((self[region.indexes] - self[region.indexes].mean())**2).sum() + 1)
