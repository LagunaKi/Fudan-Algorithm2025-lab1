def find_repeats(reference: str, query: str):
    import math

    # 构建查询序列的后缀自动机
    class SuffixAutomaton:
        def __init__(self, s: str):
            self.states = [{'link': -1, 'len': 0, 'next': {}}]
            self.last = 0
            self.count = [0]  # 每个状态的子串出现次数
            for ch in s:
                self._extend(ch)
            # 计算每个状态的endpos计数
            order = sorted(range(len(self.states)), key=lambda i: self.states[i]['len'], reverse=True)
            for v in order:
                link = self.states[v]['link']
                if link != -1:
                    self.count[link] += self.count[v]

        def _extend(self, ch: str):
            states, count = self.states, self.count
            cur = len(states)
            states.append({'link': 0, 'len': states[self.last]['len'] + 1, 'next': {}})
            count.append(1)
            p = self.last
            self.last = cur
            while p != -1 and ch not in states[p]['next']:
                states[p]['next'][ch] = cur
                p = states[p]['link']
            if p == -1:
                states[cur]['link'] = 0
            else:
                q = states[p]['next'][ch]
                if states[p]['len'] + 1 == states[q]['len']:
                    states[cur]['link'] = q
                else:
                    # 克隆状态
                    clone = len(states)
                    states.append({
                        'link': states[q]['link'],
                        'len': states[p]['len'] + 1,
                        'next': states[q]['next'].copy()
                    })
                    count.append(0)
                    while p != -1 and states[p]['next'][ch] == q:
                        states[p]['next'][ch] = clone
                        p = states[p]['link']
                    states[q]['link'] = states[cur]['link'] = clone

    # 1. 在query中构建后缀自动机并查找最长重复子串
    automaton = SuffixAutomaton(query)
    states = automaton.states
    counts = automaton.count
    max_len = 0
    max_state = 0
    for i, st in enumerate(states):
        if counts[i] >= 2 and st['len'] > max_len:
            max_len = st['len']
            max_state = i

    if max_len == 0:
        return []  # 无重复子序列

    # 还原最长重复子串
    substr = ""
    state_to_prev = {}
    for p, st in enumerate(states):
        for ch, q in st['next'].items():
            state_to_prev[q] = (p, ch)

    cur = max_state
    while states[cur]['len'] > 0:
        prev_state, ch = state_to_prev[cur]
        substr = ch + substr
        cur = prev_state

    substr = substr[-max_len:]  # 取最长匹配子串
    forward_count = counts[max_state]

    # 2. 在reference中寻找匹配（Smith-Waterman局部比对）
    def smith_waterman(ref: str, qry: str):
        n, m = len(ref), len(qry)
        dp = [[0] * (m+1) for _ in range(n+1)]
        max_score = 0
        max_pos = (0, 0)
        for i in range(1, n+1):
            for j in range(1, m+1):
                match = 2 if ref[i-1] == qry[j-1] else -1
                dp[i][j] = max(
                    0,
                    dp[i-1][j-1] + match,
                    dp[i-1][j] - 1,
                    dp[i][j-1] - 1
                )
                if dp[i][j] > max_score:
                    max_score = dp[i][j]
                    max_pos = (i, j)
        i, j = max_pos
        while i > 0 and j > 0 and dp[i][j] > 0:
            if dp[i][j] == dp[i-1][j-1] + (2 if ref[i-1] == qry[j-1] else -1):
                i -= 1; j -= 1
            elif dp[i][j] == dp[i-1][j] - 1:
                i -= 1
            else:
                j -= 1
        start = i
        length = max_pos[0] - start
        return start, length, max_score

    ref_start, substr_len, score = smith_waterman(reference, substr)
    if score == 0:
        pos = reference.find(substr)
        if pos != -1:
            ref_start = pos
            substr_len = len(substr)

    # 3. 检测反向互补序列
    comp_table = str.maketrans("ACGT", "TGCA")
    ref_rev = reference.translate(comp_table)[::-1]
    rev_flag = False
    rev_start = -1
    pos = ref_rev.find(substr)
    if pos != -1:
        rev_flag = True
        rev_start = len(reference) - pos - len(substr)

    # 4. 统计该片段在 query 中的重复次数，并确定是否发生逆转
    total_count = forward_count
    if rev_flag:
        total_count += 1

    # 5. 输出结果
    result = []
    start_position = ref_start if ref_start != -1 else rev_start
    result.append((start_position, substr_len, total_count, "是" if rev_flag else "否"))
    return result

# 测试函数（可以替换为实际的reference和query序列）
"""
ref = '''CTGCAACGTTCGTGGTTCATGTTTGAGCGATAGGCCGAAACTAACCGTGCATGCAACGTTAGTGGATCATTGTGGAACTATAGACTCAAACTAAGC
GAGCTTGCAACGTTAGTGGACCCTTTTTGAGCTATAGACGAAAACGGACCGAGGCTGCAAGGTTAGTGGATCATTTTTCAGTTTTAGACACAAACAAACCGAGC
CATCAACGTTAGTCGATCATTTTTGTGCTATTGACCATATCTCAGCGAGCCTGCAACGTGAGTGGATCATTCTTGAGCTCTGGACCAAATCTAACCGTGCCAGCA
ACGCTAGTGGATAATTTTGTTGCTATAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTTACCATCGG
ACCTCCACGAATCTGAAAAGTTTTAATTTCCGAGCGATACTTACGACCGGACCTCCACGAATCAGAAAGGGTTCACTATCCGCTCGATACATACGATCGGACCTCC
ACGACTCTGTAAGGTTTCAAAATCCGCACGATAGTTACGACCGTACCTCTACGAATCTATAAGGTTTCAATTTCCGCTGGATCCTTACGATCGGACCTCCTCGAAT
CTGCAAGGTTTCAATATCCGCTCAATGGTTACGGACGGACCTCCACGCATCTTAAAGGTTAAAATAGGCGCTCGGTACTTACGATCGGACCTCTCCGAATCTCAAA
GGTTTCAATATCCGCTTGATACTTACGATCGCAACACCACGGATCTGAAAGGTTTCAATATCCACTCTATA'''

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
"""
ref = "CTGCAACGTTCGTGGTTCATGTTTGAGCGATAGGCCGAAACTAACCGTGCATGCAACGTTAGTGGATCATTGTGGAACTATAGACTCAAACTAAGCGAGCTTGCAACGTTAGTGGACCCTTTTTGAGCTATAGACGAAAACGGACCGAGGCTGCAAGGTTAGTGGATCATTTTTCAGTTTTAGACACAAACAAACCGAGCCATCAACGTTAGTCGATCATTTTTGTGCTATTGACCATATCTCAGCGAGCCTGCAACGTGAGTGGATCATTCTTGAGCTCTGGACCAAATCTAACCGTGCCAGCAACGCTAGTGGATAATTTTGTTGCTATAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCGCTCGCTTAGCTATGGTCTATGGCGCAAAAATGATGCACTAACGAGGCAGTCTCGATTAGTGTTGGTCTATAGCAACAAAATTATCCACTAGCGTTGCTGGCTCGCTTAGCTATGGTCTATGGCGCAAAAATGATGCACTAACGAGGCAGTCTCGATTAGTGTTGGTCTATAGCAACAAAATTATCCACTAGCGTTGCTGCTTACCATCGGACCTCCACGAAT"
qry = ref  # 使用完整的参考序列作为查询序列


print(find_repeats(ref, qry)) # [(参考序列起始位置, 重复序列长度, 重复次数, 是否发生逆转)]

