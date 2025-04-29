# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 15:17:03 2018

Copyright (c) 2018 nelson<at>nanosurf.com

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""


import re
import numpy as np
import h5py
import pandas as pd
import sys, os
import time


class read:
    def __init__(self, filename=None, **kwargs):
        self.filename = filename
        self.data = None
        self.param = None

        if not filename:
            print("\nFile required!\n")
            return

        _, ext = os.path.splitext(filename)
        if ext == ".nid":
            nid = nid_read(filename=filename, **kwargs)
            self.data = nid.data
            self.param = nid.param

        elif ext == ".nhf":
            nhf = nhf_read(filename)

            self.data = nhf.data
            self.param = nhf.param

        else:
            print("Not a valid filetype")
            return


class nid_read:
    def __init__(self, filename=None, dataframe=True, verbose=False):
        self.filename = filename
        self.data = None
        self.param = None

        if not filename:
            print("\nFile required!\n")
            return

        self.__dataframe = dataframe
        self.verbose = verbose
        # names of the channels that can be opened.  Can be added to.
        self.data_names = [
            "Spec forward",
            "Spec backward",
            "Spec fwd pause",
            "Spec bwd pause",
            "Scan forward",
            "Scan backward",
            "2nd scan forward",
            "2nd scan backward",
            "Frequency sweep",
            "Frequency sweep SHO",
            "Spectrum FFT",
            "Spectrum Fit",
        ]
        self.data_types = [1, 1, 1, 1, 0, 0, 0, 0, 2, 2, 2, 2]

        self.read()

    def read(self):
        t = time.time()

        encode = "ANSI"
        # encode = 'UTF-8'
        encode = "ISO-8859-1"

        # isolates header from data file
        if self.verbose:
            print("Reading Header")

        fid = open(self.filename, "rb")
        data_in = fid.read().split(b"\n\r\n\r")
        fid.close()

        head = np.array(data_in[0].split(b"\r\n"))
        boo = head == b""
        indices = np.nonzero(boo[1:] != boo[:-1])[0] + 1
        head = np.split(head, indices)

        head = [h for h in head if len(h) > 1]

        header = {}
        for h in head:
            title = h[0].decode(encode)[1:-1]
            data = h[1:]
            values = {}
            for d in data:
                val = d.decode(encode).split("=")
                values[val[0]] = val[1].rstrip()
            header[title] = values

        param = {}

        def from_array(key, param):
            return float(header[key][param].split(",")[4])

        def get_values(key, param, df=False):
            string = header[key][param]

            ty, val, unit = re.search(
                r"(\S)\[(.*)\]\*\[(.*)\]", string
            ).groups()
            if ty == "D":
                val = float(val)
            if ty == "B":
                val = bool(val)
            if ty == "L":
                val = int(val)
            if ty == "V":
                val = [float(v) for v in val.split(",")]
                unit = [str(u) for u in unit.split(",")]
            out = {"Value": val, "Unit": unit}
            if df:
                if len(out) == 1:
                    return pd.DataFrame(out)
                else:
                    return pd.DataFrame(out)
            else:
                return out

        groupNames = []
        dataSet = []

        ds = "DataSet"
        grCnt = int(header[ds]["GroupCount"])
        for gr in range(grCnt):
            groupNames.append(header[ds]["Gr" + str(gr) + "-Name"])
            gc = int(header[ds]["Gr" + str(gr) + "-Count"])

            for g in range(gc):
                key = "Gr" + str(gr) + "-Ch" + str(g)
                if key in header[ds]:
                    dataSet.append(header[ds][key])

        terms = [
            "Frame",
            "Points",
            "Lines",
            "SaveBits",
            "Dim0Min",
            "Dim0Range",
            "Dim0Unit",
            "Dim1Min",
            "Dim1Range",
            "Dim1Unit",
            "Dim2Min",
            "Dim2Range",
            "Dim2Unit",
            "Dim2Name",
        ]

        types = [0, 1, 1, 1, 2, 2, 0, 2, 2, 0, 2, 2, 0, 0]

        # 0 = string, 1 = integer, 2 = float
        typ = [lambda x: str(x), lambda x: int(x), lambda x: float(x)]

        for term, ty in zip(terms, types):
            val = []
            for d in dataSet:
                h = header[d]
                if term in h:
                    val.append(typ[ty](h[term]))

            param[term] = val

        # for spectroscopy with maps
        dataPoints = []
        for d in dataSet:
            h = header[d]
            linePoints = []
            if "LineDim0Min" in h:
                for key, value in h.items():
                    if re.match(r"LineDim\d*Points", key):
                        linePoints.append(int(value))
            dataPoints.append(np.array(linePoints))
        param["LinePoints"] = dataPoints



        mapTable = {}
        key = "DataSet\\SpecInfos"
        if key in header:
            n = int(header[key]["SubSectionCount"])
            secNames = [header[key]["SubSection" + str(i)] for i in range(n)]

            specMode = header[key + "\\" + secNames[0]]["SpecMode"]

            count = int(header[key + "\\" + secNames[1]]["Count"])

            mapTable = [
                header[key + "\\" + secNames[1]][specMode[:3] + str(i)].split(
                    ";"
                )
                for i in range(count)
            ]
            mapTable = np.array(mapTable).astype(float)[0]

        specParam = {}
        modFw, modBw, pauseFw, pauseBw, rangeFw, rangeBw, stopFw = (
            [] for i in range(7)
        )
        key = "DataSet\\Parameters\\Spectroscopy"
        if key in header:
            modFw = get_values(key, "FwModTime")
            modBw = get_values(key, "BwModTime")
            pauseFw = get_values(key, "FwPauseTime")
            pauseBw = get_values(key, "BwPauseTime")
            rangeFw = get_values(key, "FwModRange")
            rangeBw = get_values(key, "BwModRange")
            stopFw = get_values(key, "FwModStopValue")

        scanParam = {}
        scanRange, scanTime, scanOffset, scanRotation = ([] for i in range(4))
        key = "DataSet\\Parameters\\Imaging"
        if key in header:
            scanRange = get_values(key, "ScanRange")
            scanTime = get_values(key, "ScanTime")
            scanOffset = get_values(key, "ScanOffset")
            scanRotation = get_values(key, "ScanRotation")

        thermal = {}
        key = "DataSet-Info"
        if "-- Thermal Tuning --" in header[key]:
            if "k" in header[key]["Q Factor:"]:
                Q = float(header[key]["Q Factor:"][:-1]) * 1e3
            else:
                Q = float(header[key]["Q Factor:"])
            Freq = header[key]["Frequency:"]
            SpringK = header[key]["Spring Constant:"]
            Peak = header[key]["Peak Value:"]

            gr = re.search(r"(\d+\.\d+)(\S*)", Freq).groups()
            Freq = {"Value": gr[0], "Unit": gr[1]}
            try:
                gr = re.search(r"(\d+\.\d+)\s(\S*)", SpringK).groups()
            except:
                gr = re.search(r"(\d)\s(\S*)", SpringK).groups()
            SpringK = {"Value": gr[0], "Unit": gr[1]}

            try:
                gr = re.search(r"(\d+\.\d+)(\S*)", Peak).groups()
            except:
                gr = re.search(r"(\d)(\S*)", Peak).groups()
            Peak = {"Value": gr[0], "Unit": gr[1]}

            thermal = {
                "Frequency": Freq,
                "Q Factor": Q,
                "Spring Constant": SpringK,
                "Peak": Peak,
            }

        required_size = int(
            sum(
                [
                    a * b * c / 8
                    for a, b, c in zip(
                        param["Points"], param["Lines"], param["SaveBits"]
                    )
                ]
            )
        )

        q = float(2 ** param["SaveBits"][0])
        z0 = float(q / 2)

        if param["SaveBits"][0] == 16:
            dt = np.int16
        elif param["SaveBits"][0] == 32:
            dt = np.int32

        if self.verbose:
            print("Reading Data")

        fid = open(self.filename, "rb")
        fid.seek(-required_size, 2)
        data_in = np.fromfile(fid, count=required_size, dtype=dt).astype(float)
        data_in = (data_in + z0) / q
        fid.close()

        data_in = np.split(
            data_in,
            np.cumsum(
                [a * b for a, b in zip(param["Points"], param["Lines"])]
            ),
        )
        data_in.pop()

        data = []
        # reshape and rescale data
        for datain, pts, lns, zran, zmin in zip(
            data_in,
            param["Points"],
            param["Lines"],
            param["Dim2Range"],
            param["Dim2Min"],
        ):

            data.append(datain.reshape((lns, pts)).__mul__(zran).__add__(zmin))

        data_crop = []
        try:
            for numP, dSet in zip(param["LinePoints"], data):
                if numP.any():
                    temp = []
                    for n, d in zip(numP, dSet):
                        temp.append(d[:n])
                    data_crop.append(temp)
                else:
                    data_crop.append(dSet)
            data = np.asarray(data_crop, dtype=object)
        except ValueError:
            data_crop = []
            for numP, dSet in zip(param["LinePoints"], data):
                temp = []
                # Always iterate over each row in dSet
                for i, d in enumerate(dSet):
                    # If numP has a nonzero value for this row, crop it; otherwise, use the full row.
                    if numP.any():
                        # Use the cropped version for this row
                        n = numP[i]
                        temp.append(d[:n])
                    else:
                        # Otherwise, add the entire row
                        temp.append(d)
                data_crop.append(temp)

            # Instead of forcing it into a NumPy array, keep it as a list of lists.
            #data = data_crop
            data_obj = np.empty(len(data_crop), dtype=object)
            for i, item in enumerate(data_crop):
                data_obj[i] = np.array(item, dtype=object)
            data = data_obj

        # this is the index of the channels that will be output.
        idx = [
            idx
            for idx, frame in enumerate(param["Frame"])
            if frame in self.data_names
        ]

        # only take data from approved list in name
        data = data[idx]
        frames = np.array(param["Frame"])[idx]
        channel = np.array(param["Dim2Name"])[idx]

        frame_chan = [i + ";" + j for i, j in zip(frames, channel)]

        seen = set()
        new_channel = []
        for x in frame_chan:
            if x not in seen:
                new_channel.append(x)
                seen.add(x)
            else:
                new_channel.append(x + "2")
        channel = [x.split(";")[1] for x in new_channel]

        out = {}
        for frame, chan, dat in zip(frames, channel, data):
            if frame not in out:
                out[frame] = {}
            if chan not in out[frame]:
                out[frame][chan] = {}
            if isinstance(dat, list):
                out[frame][chan] = dat
            else:
                out[frame][chan] = dat.astype(float)

        image = {}
        spec = {}
        spectrum = {}

        for name, value in out.items():
            dt = self.data_types[self.data_names.index(name)]
            if dt == 0:
                if name == "Scan forward":
                    image["Forward"] = value
                elif name == "Scan backward":
                    image["Backward"] = value
                elif name == "2nd scan forward":
                    image["2nd Forward"] = value
                elif name == "2nd scan backward":
                    image["2nd Backward"] = value
                else:
                    pass
            if dt == 1:
                if name == "Spec forward":
                    spec["Forward"] = value
                elif name == "Spec backward":
                    spec["Backward"] = value
                elif name == "Spec fwd pause":
                    spec["Pause Forward"] = value
                elif name == "Spec bwd pause":
                    spec["Pause Backward"] = value
                else:
                    pass
            if dt == 2:
                if name == "Spectrum FFT":
                    spectrum["FFT"] = value
                elif name == "Spectrum Fit":
                    spectrum["Fit"] = value
                elif name == "Frequency sweep":
                    spectrum["Sweep"] = value
                elif name == "Frequency sweep SHO":
                    spectrum["SHO"] = value
                else:
                    pass

        x = dict(
            zip(
                ["min", "range", "units"],
                [param["Dim0Min"], param["Dim0Range"], param["Dim0Unit"]],
            )
        )
        y = dict(
            zip(
                ["min", "range", "units"],
                [param["Dim1Min"], param["Dim1Range"], param["Dim1Unit"]],
            )
        )
        z = dict(
            zip(
                ["min", "range", "units"],
                [param["Dim2Min"], param["Dim2Range"], param["Dim2Unit"]],
            )
        )
        scanParam = dict(
            zip(
                ["range", "time/line", "offset", "rotation"],
                [scanRange, scanTime, scanOffset, scanRotation],
            )
        )
        specParam = dict(
            zip(
                [
                    "modFw",
                    "modBw",
                    "pauseFw",
                    "pauseBw",
                    "rangeFw",
                    "rangeBw",
                    "stop-by-value",
                    "maps",
                ],
                [
                    modFw,
                    modBw,
                    pauseFw,
                    pauseBw,
                    rangeFw,
                    rangeBw,
                    stopFw,
                    mapTable,
                ],
            )
        )

        parameters = dict(
            zip(
                ["X", "Y", "Z", "Scan", "Spec", "Tune", "HeaderDump"],
                [ x, y, z, scanParam, specParam, thermal, header],
            )
        )

        dataout = [image, spec, spectrum]
        datanames = ["Image", "Spec", "Sweep"]

        valid_data = [val for num, val in enumerate(dataout) if val]
        valid_names = [
            datanames[num] for num, val in enumerate(dataout) if val
        ]

        dataout = dict(zip(valid_names, valid_data))
        t_end = time.time() - t

        if self.verbose:
            print("Elapsed Time: %3.2f sec\n" % t_end)

        if self.__dataframe:
            # import pandas as pd
            self.data = toPandaDF(dataout).unstack(level=0)
            self.data.dropna(inplace=True)
            # self.specunits = pd.DataFrame(specUnits)
            self.param = toPandaSeries(parameters)
        else:
            self.data = dataout
            self.param = parameters


class nhf_read:
    def __init__(self, filename=None):
        self.filename = filename
        self.data = None
        self.param = None
        self.read()

    def read(self):
        # measurement
        fopen = h5py.File(self.filename, "r")

        group_name = list(fopen.keys())
        # forward/backward
        segment = fopen[group_name[0]]
        seg_names = list(segment.keys())

        s_prop = segment.attrs
        prop_s = {name: s_prop[name] for name in list(s_prop)}

        DAC_max = 2**31 - 1

        # DAC_max: dict = {0: [-(2.0**31), 2.0**31, 'dt_int32'],
        #                  1: [-(2.0**8), 2.0**8, 'dt_uint8'],
        #                  2: [-(2.0**16), 2.0**16, 'dt_uint16'],
        #                  3: [-(2.0**32), 2.0**32, 'dt_uint32'],
        #                  4: [-(2.0**64), 2.0**64, 'dt_uint64'],
        #                  5: [-(2.0**15), 2.0**15, 'dt_int16'],
        #                  6: [-(2.0**64), 2.0**64, 'dt_int64'],
        #                  7: [0.0, 0.0, 'dt_double']
        #                  }

        data_out = {}
        prop = []

        for s in seg_names:
            meas = segment[s]

            if fopen.attrs["nsf_file_version_major"] == 1:
                direction = meas.attrs["name"]
                Nx = prop_s["image_points_per_line"]
                Ny = prop_s["image_number_of_lines"]
            else:
                direction = meas.attrs["segment_name"]
                Nx, Ny = meas.attrs["rect_axis_size"]

            data_out[direction] = {}

            for m in list(meas.keys()):
                data = meas[m]
                d_prop = data.attrs

                if fopen.attrs["nsf_file_version_major"] == 1:
                    prop.append({name: d_prop[name] for name in list(d_prop)})
                    cal_max = d_prop["base_calibration_max"]
                    cal_min = d_prop["base_calibration_min"]
                    channel = d_prop["name"]

                    scale = (cal_max - cal_min) / DAC_max / 2

                    data_out[direction][channel] = (
                        np.reshape(data[:], [Nx, Ny]) * scale
                    )

                elif fopen.attrs["nsf_file_version_major"] == 2:
                    if "signal_name" in list(d_prop):
                        prop.append(
                            {name: d_prop[name] for name in list(d_prop)}
                        )
                        cal_max = d_prop["signal_calibration_max"]
                        cal_min = d_prop["signal_calibration_min"]
                        channel = d_prop["signal_name"]

                        scale = (cal_max - cal_min) / DAC_max / 2

                        data_out[direction][channel] = (
                            np.reshape(data[:], [Nx, Ny]) * scale
                        )

        _data_out = {}
        _data_out["Image"] = data_out

        self.data = toPandaDF(_data_out).unstack(level=0)
        self.param = pd.Series(prop_s)
        fopen.close()


def toPandaDF(user_dict):
    df = pd.DataFrame.from_dict(
        {
            (i, j): user_dict[i][j]
            for i in user_dict.keys()
            for j in user_dict[i].keys()
        },
        orient="index",
    ).transpose()
    return df


def toPandaSeries(user_dict):
    df = pd.Series(
        {
            (i, j): user_dict[i][j]
            for i in user_dict.keys()
            for j in user_dict[i].keys()
        }
    )
    return df


def checkFileType(filename, cla):
    import sys, os

    switch = {"nid": "nhf_read", "nhf": "nid_read"}

    _, ext = os.path.splitext(filename)
    if ext[-3:] != cla:
        print(f'\nWrong file type detected, use "%s" instead' % switch[cla])
        sys.exit()
    elif ext[-3:] == cla:
        return filename
    else:
        print("\nNot a valid filetype")
        exit()


if __name__ == "__main__":
    # filename='..//example//Plane_Removal//circles.nid'
    filename = "atomic_lattice.nid"
    afm = read(filename, verbose=True)
    print(afm.param)
