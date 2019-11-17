
def bin2dec(bits, length):
    if type(bits) == str:
        bits = list(map(int, list(bits)))
    bits = bits[::-1]
    result = 0
    if bits[length - 1] == 0:
        for i in range(length):
            result += bits[i] * (2 ** i)
    else:
        for i in range(length - 1):
            result += bits[i] * (2 ** i)
        result = result - 2 ** (length - 1)

    return result

def find(seq):
  """Return first item in sequence where f(item) == True."""
  f = lambda number: number != None
  for item in seq:
    if f(item):
      return item

def matrix_sum(oper, M, N):
    result = 0

    for i in range(M):
        oper[i] = [find(oper[i]) if oper[i][j] == None else oper[i][j] for j in range(len(oper[i]))]
        oper[i] += [0] * (N - len(oper[i]))

    for i in range(M):
        result += bin2dec(oper[i], N)

    return result