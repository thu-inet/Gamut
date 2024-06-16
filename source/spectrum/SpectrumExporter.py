import xml.etree.ElementTree as ET
import pickle
import pandas as pd

from time import strftime, localtime
from pathlib import Path
from xml.dom.minidom import parse as parse_xml


class SpectrumExporter:

    @staticmethod
    def export_to_pickle(spectrum, filename: str | Path) -> None:
        """
        Export to binary pickle file.
        
        :param spectrum: Spectrum object to be exported.
        :param filename: Path to the file where the spectrum will be saved.
        """
        with open(filename, 'wb') as fopen:
            pickle.dump(spectrum, fopen)

    @staticmethod
    def export_to_GammaVision(spectrum, filename: str | Path) -> None:
        """
        Export to Ortec plain text file.

        :param spectrum: Spectrum object to be exported.
        :param filename: Path to the file where the spectrum will be saved.
        """
        with open(filename, 'w') as fopen:
            fopen.write("$SPEC_ID:\n\n")
            fopen.write("$SPEC_REM:\n\n")
            fopen.write("$DATE_MEA:\n" + strftime(r"%m/%d/%Y %H:%M:%S", localtime()) + "\n")
            fopen.write("$MEAS_TIME:\n" + "1000 1000\n")
            fopen.write(f"$DATA:\n0 {spectrum.length-1}\n")
            fopen.write("\n".join([f"{round(c):>8d}" for c in spectrum.counts]) + "\n")
            if len(spectrum.regions) != 0:
                fopen.write(f"$ROI:\n{len(spectrum.regions)}\n" + "\n".join([f"{region.left} {region.right}" for region in spectrum.regions]))
            fopen.write("$PRESETS:\nNone\n0\n0\n")
            fopen.write(f"$ENER_FIT:\n{spectrum.ergcal.params[0]:8.6f} {spectrum.ergcal.params[1]:8.6f}\n")

    @staticmethod
    def export_to_pandas(spectrum) -> pd.DataFrame:
        """
        Export to excel file.
        
        :param spectrum: Spectrum object to be exported.
        
        :return: DataFrame with the ROIs data.
        """
        results = []
        for region in spectrum.regions:
            for peak in region.peaks:

                val = {'left': region.left, 'right': region.right,
                        'energy': spectrum.ergcal(peak['location'])}
                attrs = ['location', 'height', 'stderr', 'area', 'sig_area2', 'fitness']
                val.update(dict([(attr, peak[attr]) for attr in attrs if attr in peak.keys()]))
                if 'stderr' in peak.keys():
                    val.update({'fitness': spectrum.estimate_fitness(region)})  

                results.append(val)
        results = pd.DataFrame(results)
        return results
    
    @staticmethod
    def export_to_xml(spectrum, filepath: str | Path,
                    data: dict[str, bool] = {'ergcal': True, 'ergcal.verbose': False,
                                            'FWHMcal': True, 'FWHMcal.verbose': False,
                                            'regions': True, 'regions.verbose': False,
                                            'regions.peaks': True, 'regions.peaks.verbose': False}) -> None:
        """
        Export the report to a xml file.
        
        :param spectrum: Spectrum object to be exported.
        :param filename: Path to the file where the spectrum will be saved.
        :param data: Dictionary indicating which data will be included in the report.
        """
        spec = ET.Element('Spectrum')
        spec.attrib['name'] = spectrum.label
        spec.attrib['length'] = str(spectrum.length)

        if data['ergcal']:
            ergcal_el = ET.SubElement(spec, 'ergcal')
            ergcal_el.attrib['method'] = spectrum.ergcal.method
            ergcal_el.attrib['params'] = " ".join([f"{val:.4f}" for val in spectrum.ergcal.params])
            if data['ergcal.verbose'] and hasattr(spectrum.ergcal, 'data'):
                ET.SubElement(spec, 'data_indexes').text = " ".join([f"{val:.2f}" for val in spectrum.ergcal.data_indexes])
                ET.SubElement(spec, 'data_values').text = " ".join([f"{val:.2f}" for val in spectrum.ergcal.data_values])
                ergcal_el.attrib['fitness'] = f"{spectrum.ergcal.fitness:.2f}"

        if data['FWHMcal']:
            FWHMcal_el = ET.SubElement(spec, 'FWHMcal')
            FWHMcal_el.attrib['method'] = spectrum.FWHMcal.method
            FWHMcal_el.attrib['params'] = " ".join([f"{val:.4f}" for val in spectrum.FWHMcal.params])
            if data['FWHMcal.verbose'] and hasattr(spectrum.FWHMcal, 'data'):
                ET.SubElement(spec, 'data_ind').text = " ".join([f"{val:.2f}" for val in spectrum.FWHMcal.data[0, :]])
                ET.SubElement(spec, 'data_val').text = " ".join([f"{val:.2f}" for val in spectrum.FWHMcal.data[1, :]])
                FWHMcal_el.attrib['fitness'] = f"{spectrum.FWHMcal.fitness:.2f}"

        if data['regions']:
            regions_el = ET.SubElement(spec, 'regions')
            for region in spectrum.regions:
                region_el = ET.SubElement(regions_el, 'region')
                region_el.attrib['N_peaks'] = str(len(region.peaks))
                region_el.attrib['ind_left'] = str(region.left)
                region_el.attrib['ind_right'] = str(region.right)

                if data['regions.verbose']:
                    region_el.attrib['erg_left'] = f"{spectrum.ergcal(region.left):.2f}"
                    region_el.attrib['erg_right'] = f"{spectrum.ergcal(region.right):.2f}"
                    try:
                        region_el.attrib['fitness'] = f"{spectrum.estimate_fitness(region):.2f}"
                        region_el.attrib['slope'] = f"{region.slope:.2f}"
                        region_el.attrib['offset'] = f"{region.offset:.2f}"
                    except:
                        pass
                
                if data['regions.peaks']:
                    for peak in region.peaks:
                        peak_el = ET.SubElement(region_el, 'peak')
                        peak_el.attrib['energy'] = f"{spectrum.ergcal(peak['location']):.2f}"
                        for key, val in peak.items():
                            peak_el.attrib[key] = str(val) if isinstance(val, int) else f"{val:.2f}"

        ET.SubElement(spec, 'counts').text = " ".join([f"{val:.2f}" for val in spectrum.counts])

        tree = ET.ElementTree(spec)
        tree.write(filepath, encoding='utf-8', xml_declaration=True)
        dom = parse_xml(str(filepath))
        pdom = dom.toprettyxml(indent='\t', newl='\n')
        dom.unlink()
        with open(filepath, 'w') as fileopen:
            fileopen.write(pdom)