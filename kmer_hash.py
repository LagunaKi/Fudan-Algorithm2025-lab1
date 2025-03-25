from mysort import my_sort

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[i] for i in reversed(seq))

# 时间复杂度O(k)
def build_kmer_hash_table(S, k):
    kmer_counts = {}
    n = len(S)
    # 因为是reference的哈希表，所以只考虑每个kmer第一次出现的位置
    for i in range(n - k + 1):
        kmer = S[i:i + k]
        rev_comp = reverse_complement(kmer)
        if kmer not in kmer_counts:
            kmer_counts[kmer] = i + k
        if rev_comp not in kmer_counts:
            kmer_counts[rev_comp] = i + k
    return kmer_counts

def find_repeats_kmer(ref, qry, max_len=100, min_len=50):
    repeats = []
    m = len(qry)
    covered = [False] * m  # 标记检测到的重复区域

    for k in range(max_len, min_len - 1, -1): # 先查找长的重复
        if k > len(ref) or k > len(qry):
            continue

        hash_table = build_kmer_hash_table(ref, k) # 构建reference的kmer哈希表

        # 遍历query序列
        j = 0
        while j <= len(qry) - k:
            if covered[j]:
                j += 1
                continue
            
            current_kmer = qry[j:j + k]

            if current_kmer in hash_table:
                ref_pos = hash_table[current_kmer]
                is_reverse = (current_kmer != ref[ref_pos - k:ref_pos])
                
                # 寻找连续重复，统计重复次数
                repeat_count = 1
                next_j = j + k
                while next_j <= len(qry) - k:
                    next_kmer = qry[next_j:next_j + k]
                    if next_kmer == current_kmer:
                        repeat_count += 1
                        next_j += k
                    else:
                        break

                if repeat_count > 1 :
                    repeats.append((ref_pos, k, repeat_count, is_reverse, j))
                    # 标记已检测的区域
                    for pos in range(j, next_j):
                        if pos < m:
                            covered[pos] = True
                # 跳到未检测的区域
                j = next_j
            else:
                j += 1


    # 处理连续重复序列（第一个重复的count减一）
    my_sort(repeats, key=lambda x: x[4]) # 按j排序，自定义快排
    adjusted_repeats = []
    i = 0
    while i < len(repeats):
        current = repeats[i]
        if i < len(repeats)-1:
            next_rep = repeats[i+1] # 下一个重复
            current_end = current[4] + current[1] * current[2] # j + k * count, 结束位置
            # 判断是否连续且类型不同
            if (next_rep[4] == current_end and (current[0], current[1], current[3]) != (next_rep[0], next_rep[1], next_rep[3])):
                adjusted = (current[0], current[1], current[2]-1, current[3], current[4])
                adjusted_repeats.append(adjusted)
                adjusted_repeats.append(next_rep)
                i += 2
                continue
        adjusted_repeats.append(current)
        i += 1
    repeats = adjusted_repeats
    
    # 去重处理
    seen = set()
    unique_repeats = []
    for pos, length, count, rev, _ in repeats:
        key = (pos, length, rev)
        if key not in seen:
            seen.add(key)
            unique_repeats.append((pos, length, count, rev))
    
    return unique_repeats

def read_data(file_path):
    with open(file_path, 'r') as file:
        return file.read().replace('\n', '')

reference = read_data('reference.txt')
query = read_data('query.txt')
repeats = find_repeats_kmer(reference, query)

for idx, (pos, length, count, is_reverse) in enumerate(repeats, 1):
    print(f"重复{idx}: 位置={pos}, 长度={length}, 重复次数={count}, 反向={is_reverse}")