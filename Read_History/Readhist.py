# This reads a .hst midas history file and makes a plot. Useful to perform analysis on history data
# This is a template and can be further complicated reading multiple variables, multiple files and complicate plots
import struct
import argparse
import matplotlib.pyplot as plt
import numpy as np


RT_DEF  = 0x46445348
RT_DATA = 0x41445348

tid_size = {
    1:1,2:1,3:1,
    4:2,5:2,
    6:4,7:4,8:4,9:4,
    10:8,11:1
}

def main():
    parser = argparse.ArgumentParser(
    description="Read a history midas file."
    )
    parser.add_argument("--file", required=True,
                    help="Filename, es: ~/flash-data/data/LNF/260202.hst")
    parser.add_argument("--event-id", required=True,
                    help="Event id of the equipment to read, es: 3")
    parser.add_argument("--var", default='VGEM3 2 CurrentDet',
                    help="name of the variable to plot, es: VGEM3 2 CurrentDet")                  
    args = parser.parse_args()

    defs = {}
    filename = args.file
    event_id_look = int(args.event_id)
    variable = args.var

    time = np.array([])
    y = np.array([])

    with open(filename,"rb") as f:

        while True:

            h = f.read(20)

            if len(h) < 20:
                break

            record_type,event_id,timestamp,def_offset,data_size = struct.unpack("<IIIII",h)

            if record_type == RT_DEF:

                event_name = f.read(32).split(b"\0")[0].decode()

                consumed = 32
                offset = 0
                tags = []

                while consumed < data_size:

                    raw = f.read(40)

                    name = raw[:32].split(b"\0")[0].decode()
                    typ,ndata = struct.unpack("<II",raw[32:])

                    tags.append((name,typ,ndata,offset))
                    offset += tid_size.get(typ,0)*ndata

                    consumed += 40

                defs[event_id]=tags

            elif record_type == RT_DATA:

                blob = f.read(data_size)

                if event_id != event_id_look:
                    continue

                #print(f"\nTime {timestamp}")
                time = np.append(time,timestamp)

                for name,typ,ndata,offset in defs[event_id_look]:

                    if typ == 6:
                        value = struct.unpack_from("<I",blob,offset)[0]

                    elif typ == 7:
                        value = struct.unpack_from("<i",blob,offset)[0]

                    elif typ == 9:
                        value = struct.unpack_from("<f",blob,offset)[0]

                    elif typ == 10:
                        value = struct.unpack_from("<d",blob,offset)[0]

                    else:
                        value = f"type {typ}"

                    if name == variable:
                        y = np.append(y,value)
                    #print(name,value)

            else:
                f.seek(data_size,1)
    
    #filter here time if you want
    
    fig, axes1 =  plt.subplots()
    axes1.plot(time, y)
    #axes1.semilogy(time, y)
    #axes1.set_ylim(1.786*10**8,1.788*10**8)
    axes1.set(title=variable,xlabel='UNIX t (s)',ylabel='Current (uA)')
    plt.show()


if __name__ == "__main__":
    main()
