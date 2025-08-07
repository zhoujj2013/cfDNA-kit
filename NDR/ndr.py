import pysam
import numpy as np

def ndr_signal(sf:pysam.AlignmentFile,chrid,tss_start):
    # sf = pysam.AlignmentFile(bam_file, "rb")
    binedges = range(-1000, 1000, 5)
    tss = tss_start  # TSS 位置（start）
    s = tss - 1000
    e = tss + 1000

    pos = []
    iter = sf.fetch(chrid, s, e)  # 提取该区域的配对 reads
    for x in iter:
        flaglen = abs(int(x.template_length))  # 获取插入片段长度

        if x.flag == 99 and x.next_reference_id == x.reference_id:
            if flaglen < 2000:
                xs = x.reference_start
                xe = x.reference_start + flaglen
                for p in range(xs, xe):
                    rpos = p - tss  # 相对 TSS 的距离
                    if -1000 < rpos < 1000:
                        pos.append(rpos)
        elif x.flag == 83 and x.next_reference_id == x.reference_id:
            if flaglen < 2000:
                xs = x.reference_start
                xe = x.reference_start - flaglen
                for p in range(xe, xs):
                    rpos = p - tss
                    if -1000 < rpos < 1000:
                        pos.append(rpos)

    pos.append(0)  # 避免空数组报错
    vect_pos = np.histogram(pos, binedges)  # 统计直方图（TSS ±1kb）

    ### ========= 上游（-3k ~ -1k）平均深度 =========
    utss = tss_start-2000
    us = tss_start - 3000
    ue = tss_start - 1000
    upos = []
    iter = sf.fetch(chrid, us, ue)
    for x in iter:
        flaglen = abs(int(x.template_length))
        if x.flag == 99 and x.next_reference_id == x.reference_id:
            if flaglen < 2000:
                xs = x.reference_start
                xe = x.reference_start + flaglen
                for p in range(xs, xe):
                    rpos = p - utss
                    if -1000 < rpos < 1000:
                        upos.append(rpos)
        elif x.flag == 83 and x.next_reference_id == x.reference_id:
            if flaglen < 2000:
                xs = x.reference_start
                xe = x.reference_start - flaglen
                for p in range(xe, xs):
                    rpos = p - utss
                    if -1000 < rpos < 1000:
                        upos.append(rpos)

    upos.append(0)
    vect_upos = np.histogram(upos, binedges)
    upMean = np.mean(vect_upos[0])  # 上游平均深度

    ### ========= 下游（+1k ~ +3k）平均深度 =========
    dtss = tss_start + 2000
    ds = tss_start + 1000
    de = tss_start + 3000
    dpos = []
    iter = sf.fetch(chrid, ds, de)
    for x in iter:
        flaglen = abs(int(x.template_length))
        if x.flag == 99 and x.next_reference_id == x.reference_id:
            if flaglen < 2000:
                xs = x.reference_start
                xe = x.reference_start + flaglen
                for p in range(xs, xe):
                    rpos = p - dtss
                    if -1000 < rpos < 1000:
                        dpos.append(rpos)
        elif x.flag == 83 and x.next_reference_id == x.reference_id:
            if flaglen < 2000:
                xs = x.reference_start
                xe = x.reference_start - flaglen
                for p in range(xe, xs):
                    rpos = p - dtss
                    if -1000 < rpos < 1000:
                        dpos.append(rpos)

    dpos.append(0)
    vect_dpos = np.histogram(dpos, binedges)
    downMean = np.mean(vect_dpos[0])  # 下游平均深度

    updownMean = (upMean + downMean) / 2  # 上下游平均深度
    vectt_pos = [float(depth / updownMean) for depth in vect_pos[0]]  # 归一化
    return vectt_pos