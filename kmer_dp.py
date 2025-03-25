from collections import defaultdict

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

def build_kmer_hash_table(S, k):
    kmer_counts = defaultdict(list)
    n = len(S)
    for i in range(n - k + 1):
        kmer = S[i:i + k]
        rev_comp = reverse_complement(kmer)
        kmer_counts[kmer].append(i)
        kmer_counts[rev_comp].append(i)
    # 只保留出现一次的kmer
    unique_kmers = {}
    for km, positions in kmer_counts.items():
        if len(positions) == 1:
            unique_kmers[km] = positions[0]
    return unique_kmers

def find_repeats_kmer(S_ref, S_query, max_len=100, min_len=50):
    repeats = []
    m = len(S_query)
    covered = [False] * m  # 标记已覆盖的位置

    # 从最大k开始处理，优先检测长重复
    for k in range(max_len, min_len - 1, -1):
        if k > len(S_ref) or k > len(S_query):
            continue
        hash_table = build_kmer_hash_table(S_ref, k)
        if not hash_table:
            continue
        
        j = 0
        while j <= len(S_query) - k:
            if covered[j]:
                j += 1
                continue
            
            current_kmer = S_query[j:j + k]
            if current_kmer in hash_table:
                ref_pos = hash_table[current_kmer]
                is_reverse = current_kmer != S_ref[ref_pos:ref_pos + k]
                
                # 统计连续重复次数
                repeat_count = 1
                next_j = j + k
                while next_j <= len(S_query) - k:
                    next_kmer = S_query[next_j:next_j + k]
                    if next_kmer == current_kmer:
                        repeat_count += 1
                        next_j += k
                    else:
                        break
                
                # 计算额外重复次数
                extra_count = repeat_count - 1
                if extra_count > 0:
                    repeats.append((ref_pos, k, extra_count, is_reverse))
                    # 标记覆盖区域
                    for pos in range(j, next_j):
                        if pos < m:
                            covered[pos] = True
                # 跳到未覆盖区域
                j = next_j
            else:
                j += 1
    
    # 去重处理
    seen = set()
    unique_repeats = []
    for pos, length, count, rev in repeats:
        key = (pos, length, rev)
        if key not in seen:
            seen.add(key)
            unique_repeats.append((pos, length, count, rev))
    
    return unique_repeats

# 示例用法
def read_data(file_path):
    with open(file_path, 'r') as file:
        return file.read().replace('\n', '')

reference = read_data('reference.txt')
query = read_data('query.txt')
repeats = find_repeats_kmer(reference, query, max_len=100, min_len=50)

for idx, (pos, length, count, is_reverse) in enumerate(repeats, 1):
    print(f"重复{idx}: 位置={pos}, 长度={length}, 重复次数={count}, 反向={is_reverse}")