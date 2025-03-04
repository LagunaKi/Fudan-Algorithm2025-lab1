def build_suffix_automaton(s: str):
    """构建字符串 s 的后缀自动机，返回自动机的转移、后缀链接、状态最长长度、每个状态endpos计数、以及每个状态的一个示例end位置。"""
    # 初始化
    trans = [{}]           # 状态0转移
    link = [-1]            # 后缀链接
    max_len = [0]          # 状态的最长子串长度
    occ_count = [0]        # endpos计数（初始为0，稍后计算）
    endpos_example = [-1]  # 状态的一个end位置示例
    last = 0               # last指向当前尾部状态
    # 构建自动机
    for i, ch in enumerate(s):
        cur = len(trans)    # 新状态索引
        trans.append({})    
        max_len.append(max_len[last] + 1)
        link.append(0)      # 临时设为0（root），稍后可能调整
        occ_count.append(1) # 每新增一个字符，形成一个新的end
        endpos_example.append(i)
        p = last
        # 尝试从last状态开始添加转移
        while p != -1 and ch not in trans[p]:
            trans[p][ch] = cur
            p = link[p]
        if p == -1:
            # 没有找到带ch转移的前缀，当前状态的后缀链接保持为0（root）
            link[cur] = 0
        else:
            q = trans[p][ch]
            if max_len[p] + 1 == max_len[q]:
                # 情况1：可以直接将cur的链接指向q
                link[cur] = q
            else:
                # 情况2：需要克隆状态q'
                clone = len(trans)
                trans.append(trans[q].copy())   # 复制q的转移
                max_len.append(max_len[p] + 1)
                link.append(link[q])            # 克隆的链接与q相同
                occ_count.append(0)             # 克隆状态初始不赋值计数，在后面汇总
                endpos_example.append(endpos_example[q])
                # 更新p链上所有指向q的转移为指向clone
                while p != -1 and trans[p].get(ch) == q:
                    trans[p][ch] = clone
                    p = link[p]
                # 调整q和cur的后缀链接
                link[q] = clone
                link[cur] = clone
        last = cur
    # 通过后缀链接传播occ_count（从长串状态向短串状态累加）
    order = sorted(range(len(trans)), key=lambda x: max_len[x], reverse=True)
    for state in order:
        if link[state] != -1:
            occ_count[link[state]] += occ_count[state]
            # 传播endpos示例：如果链接状态还没有示例，则继承
            if endpos_example[link[state]] == -1:
                endpos_example[link[state]] = endpos_example[state]
    return trans, link, max_len, occ_count, endpos_example

def filter_and_merge_repeats(results):
    """合并连续重复的片段，仅保留最长的片段"""
    if not results:
        return []

    # 按起始位置排序
    results.sort(key=lambda x: (x[0], -x[1]))  # 先按起始位置升序，再按长度降序

    filtered_results = []
    prev_start, prev_length, prev_count, prev_reversed = results[0]

    for i in range(1, len(results)):
        start, length, count, reversed_flag = results[i]

        # 只保留最长的不重叠片段
        if start <= prev_start + prev_length - 1:  # 发现连续或嵌套
            continue  # 跳过当前短片段
        else:
            filtered_results.append((prev_start, prev_length, prev_count, prev_reversed))
            prev_start, prev_length, prev_count, prev_reversed = start, length, count, reversed_flag

    # 别忘了最后一个片段
    filtered_results.append((prev_start, prev_length, prev_count, prev_reversed))
    return filtered_results


def find_repeats(reference: str, query: str, min_length: int = 5):
    """在reference和query中寻找精确匹配的重复片段，检测反向互补出现，输出结果表格。"""
    trans, link, max_len, occ_count, endpos_ex = build_suffix_automaton(query)
    results = []         # 收集结果 (pos, length, count, reversed_flag)
    seen = set()         # 用于避免重复报告相同片段或其反向互补

    def rev_comp(seq: str) -> str:
        comp_map = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
        return "".join(comp_map[ch] for ch in seq[::-1])

    for state, count in enumerate(occ_count):
        if state == 0 or count < 2:
            continue  
        L = max_len[state]
        if L < min_length:
            continue  
        end_idx = endpos_ex[state]
        subseq = query[end_idx - L + 1 : end_idx + 1]
        comp_subseq = rev_comp(subseq)
        if subseq != comp_subseq:
            cur = 0
            for ch in comp_subseq:
                if ch in trans[cur]:
                    cur = trans[cur][ch]
                else:
                    cur = -1
                    break
            if cur != -1:
                continue  
        pos = reference.find(subseq)
        if pos == -1:
            pos = reference.find(comp_subseq)
            if pos == -1:
                continue  
        key = tuple(sorted([subseq, comp_subseq]))
        if key in seen:
            continue
        seen.add(key)
        results.append((pos, L, count, False))

    rev_query = rev_comp(query)
    lcs_len = [0] * len(trans)
    cur_state = 0
    cur_length = 0
    for ch in rev_query:
        if ch in trans[cur_state]:
            cur_state = trans[cur_state][ch]
            cur_length += 1
        else:
            while cur_state != -1 and ch not in trans[cur_state]:
                cur_state = link[cur_state]
            if cur_state == -1:
                cur_state = 0
                cur_length = 0
            else:
                cur_length = max_len[cur_state] + 1
                cur_state = trans[cur_state][ch]
        if cur_length > lcs_len[cur_state]:
            lcs_len[cur_state] = cur_length

    for state in range(len(trans)):
        if link[state] != -1:
            lcs_len[link[state]] = max(lcs_len[link[state]], min(lcs_len[state], max_len[link[state]]))

    for state, L in enumerate(lcs_len):
        if state == 0 or L < min_length:
            continue
        end_idx = endpos_ex[state]
        subseq = query[end_idx - L + 1 : end_idx + 1]
        comp_subseq = rev_comp(subseq)
        if subseq == comp_subseq:
            continue
        key = tuple(sorted([subseq, comp_subseq]))
        if key in seen:
            continue
        occ_fwd = occ_count[state]
        cur = 0
        for ch in comp_subseq:
            if ch in trans[cur]:
                cur = trans[cur][ch]
            else:
                cur = -1
                break
        occ_rev = occ_count[cur] if cur != -1 else 0
        total_count = occ_fwd + occ_rev
        if total_count < 2:
            continue  
        pos = reference.find(subseq)
        if pos == -1:
            pos = reference.find(comp_subseq)
            if pos == -1:
                continue
        seen.add(key)
        results.append((pos, L, total_count, True))

    results_by_pos = {}
    for pos, length, count, rev_flag in results:
        if pos not in results_by_pos or length > results_by_pos[pos][1]:
            results_by_pos[pos] = (pos, length, count, rev_flag)
    final_results = sorted(results_by_pos.values(), key=lambda x: x[0])

    # **应用去重和合并策略**
    final_results = filter_and_merge_repeats(final_results)

    print(f"{'Reference_Pos':>12} {'Length':>6} {'Count':>6} {'Reversed':>8}")
    for pos, length, count, rev_flag in final_results:
        print(f"{pos:12d} {length:6d} {count:6d} {'Yes' if rev_flag else 'No':>8s}")


# 测试函数
ref = '''CTGCAACGTTCGTGGTTCATGTTTGAGCGATAGGCCGAAACTAACCGTGCATGCAACGTTAGTGGATCATTGTGGAACTATAGACTCAAACTAAGCGA
GCTTGCAACGTTAGTGGACCCTTTTTGAGCTATAGACGAAAACGGACCGAGGCTGCAAGGTTAGTGGATCATTTTTCAGTTTTAGACACAAACAAACCGAGCCATCA
ACGTTAGTCGATCATTTTTGTGCTATTGACCATATCTCAGCGAGCCTGCAACGTGAGTGGATCATTCTTGAGCTCTGGACCAAATCTAACCGTGCCAGCAACGCTAG
TGGATAATTTTGTTGCTATAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTTACCATCGGACCTCCAC
GAATCTGAAAAGTTTTAATTTCCGAGCGATACTTACGACCGGACCTCCACGAATCAGAAAGGGTTCACTATCCGCTCGATACATACGATCGGACCTCCACGACTCTG
TAAGGTTTCAAAATCCGCACGATAGTTACGACCGTACCTCTACGAATCTATAAGGTTTCAATTTCCGCTGGATCCTTACGATCGGACCTCCTCGAATCTGCAAGGTT
TCAATATCCGCTCAATGGTTACGGACGGACCTCCACGCATCTTAAAGGTTAAAATAGGCGCTCGGTACTTACGATCGGACCTCTCCGAATCTCAAAGGTTTCAATAT
CCGCTTGATACTTACGATCGCAACACCACGGATCTGAAAGGTTTCAATATCCACTCTATA
'''

qry =  '''CTGCAACGTTCGTGGTTCATGTTTGAGCGATAGGCCGAAACTAACCGTGCATGCAACGTTAGTGGATCATTGTGGAACTATAGACTCAAACTAAGCG
AGCTTGCAACGTTAGTGGACCCTTTTTGAGCTATAGACGAAAACGGACCGAGGCTGCAAGGTTAGTGGATCATTTTTCAGTTTTAGACACAAACAAACCGAGCCATC
AACGTTAGTCGATCATTTTTGTGCTATTGACCATATCTCAGCGAGCCTGCAACGTGAGTGGATCATTCTTGAGCTCTGGACCAAATCTAACCGTGCCAGCAACGCTA
GTGGATAATTTTGTTGCTATAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATC
ATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTG
CGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAG
TGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAG
CTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCGCTCGCTTAGCTATGGTCTATGGCGCAAAAATGATG
CACTAACGAGGCAGTCTCGATTAGTGTTGGTCTATAGCAACAAAATTATCCACTAGCGTTGCTGGCTCGCTTAGCTATGGTCTATGGCGCAAAAATGATGCACTAAC
GAGGCAGTCTCGATTAGTGTTGGTCTATAGCAACAAAATTATCCACTAGCGTTGCTGCTTACCATCGGACCTCCACGAATCTGAAAAGTTTTAATTTCCGAGCGATA
CTTACGACCGGACCTCCACGAATCAGAAAGGGTTCACTATCCGCTCGATACATACGATCGGACCTCCACGACTCTGTAAGGTTTCAAAATCCGCACGATAGTTACGA
CCGTACCTCTACGAATCTATAAGGTTTCAATTTCCGCTGGATCCTTACGATCGGACCTCCTCGAATCTGCAAGGTTTCAATATCCGCTCAATGGTTACGGACGGACC
TCCACGCATCTTAAAGGTTAAAATAGGCGCTCGGTACTTACGATCGGACCTCTCCGAATCTCAAAGGTTTCAATATCCGCTTGATACTTACGATCGCAACACCACGG
ATCTGAAAGGTTTCAATATCCACTCTATA
'''



ref = ref.replace('\n', '')  # 移除所有换行符
qry = qry.replace('\n', '')




find_repeats(ref, qry, min_length=15)

