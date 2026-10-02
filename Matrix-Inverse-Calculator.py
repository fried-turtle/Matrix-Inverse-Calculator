def matrixout(mx, size):
    formatted_rows = []
    for i in range(size):
        row_str = ""
        for j in range(size):
            val = mx[i][j]
            if abs(val - round(val)) < 1e-6:
                row_str += "%9d" % int(round(val))
            else:
                row_str += "%9.4f" % val
        formatted_rows.append(row_str)

    content_width = len(formatted_rows[0])

    print("┌" + " " * content_width + "┐")
    for row in formatted_rows:
        print("│" + row + "│")
    print("└" + " " * content_width + "┘")


def transpose(m):
    return [
        [m[j][i] for j in range(len(m))]
        for i in range(len(m[0]))
    ]


def get_minor(m, i, j):
    return [
        row[:j] + row[j + 1:]
        for row in (m[:i] + m[i + 1:])
    ]


def get_det(m):
    n = len(m)

    if n == 1:
        return m[0][0]

    if n == 2:
        return m[0][0] * m[1][1] - m[0][1] * m[1][0]

    determinant = 0
    for c in range(n):
        minor = get_minor(m, 0, c)
        determinant += ((-1) ** c) * m[0][c] * get_det(minor)

    return determinant


def inv_by_det(m):
    determinant = get_det(m)

    if abs(determinant) < 1e-10:
        return None, "행렬식이 0이므로 역행렬이 존재하지 않습니다."

    n = len(m)

    if n == 1:
        return [[1.0 / m[0][0]]], None

    if n == 2:
        inverse = [
            [m[1][1] / determinant, -m[0][1] / determinant],
            [-m[1][0] / determinant, m[0][0] / determinant]
        ]
        return inverse, None

    cofactor_matrix = []
    for r in range(n):
        cofactor_row = []
        for c in range(n):
            minor = get_minor(m, r, c)
            cofactor = ((-1) ** (r + c)) * get_det(minor)
            cofactor_row.append(cofactor)
        cofactor_matrix.append(cofactor_row)

    adjugate = transpose(cofactor_matrix)

    for r in range(n):
        for c in range(n):
            adjugate[r][c] /= determinant

    return adjugate, None


def inv_by_gauss(m):
    n = len(m)
    aug = []

    for i in range(n):
        identity_row = [1.0 if i == j else 0.0 for j in range(n)]
        aug.append(m[i][:] + identity_row)

    for i in range(n):
        max_row = i
        for r in range(i + 1, n):
            if abs(aug[r][i]) > abs(aug[max_row][i]):
                max_row = r

        aug[i], aug[max_row] = aug[max_row], aug[i]
        pivot = aug[i][i]

        if abs(pivot) < 1e-10:
            return None, "역행렬이 존재하지 않는 특이행렬입니다."

        for c in range(2 * n):
            aug[i][c] /= pivot

        for r in range(n):
            if r != i:
                factor = aug[r][i]
                for c in range(2 * n):
                    aug[r][c] -= factor * aug[i][c]

    inverse = [row[n:] for row in aug]
    return inverse, None


def compare_matrix(a, b, size):
    tolerance = 1e-4
    for i in range(size):
        for j in range(size):
            if abs(a[i][j] - b[i][j]) > tolerance:
                return False
    return True


def multiply_matrix(a, b, size):
    result = [[0.0 for _ in range(size)] for _ in range(size)]
    for i in range(size):
        for j in range(size):
            for k in range(size):
                result[i][j] += a[i][k] * b[k][j]
    return result


def verify_inverse(a, inverse, size):
    print("\n[추가기능] 역행렬 수학적 검증 (A × A⁻¹ = I)")
    result = multiply_matrix(a, inverse, size)

    print("\n[A × A⁻¹]")
    matrixout(result, size)

    tolerance = 1e-4
    for i in range(size):
        for j in range(size):
            expected = 1.0 if i == j else 0.0
            if abs(result[i][j] - expected) > tolerance:
                print("\n=> 검증 결과: 실패 (단위행렬 미일치)")
                return False

    print("\n=> 검증 결과: 성공 (A × A⁻¹ = I 도출)")
    return True


def input_matrix(size):
    matrix = []
    for i in range(size):
        while True:
            try:
                row = [float(x) for x in input(f"{i + 1}행: ").split()]
                if len(row) != size:
                    print(f"입력 오류: 정확히 {size}개의 값을 입력해야 합니다.")
                    continue
                matrix.append(row)
                break
            except ValueError:
                print("입력 오류: 숫자만 입력해주세요.")
    return matrix


def main():
    try:
        size = int(input("정방행렬의 차수를 입력하세요: "))
        if size <= 0:
            print("오류: 차수는 양의 정수여야 합니다.")
            return

        matrix = input_matrix(size)

        # 1. 행렬식 방식
        inv_det, err_det = inv_by_det(matrix)
        print("\n1. 행렬식으로 구한 역행렬:")
        if err_det:
            print("오류:", err_det)
        else:
            matrixout(inv_det, size)

        # 2. 가우스-조던 방식
        inv_gj, err_gj = inv_by_gauss(matrix)
        print("\n2. 가우스-조던 소거법으로 구한 역행렬:")
        if err_gj:
            print("오류:", err_gj)
        else:
            matrixout(inv_gj, size)

        # 3. 결과 비교
        print("\n3. 결과 비교:")
        if err_det and err_gj:
            print("두 방법 모두 역행렬을 구할 수 없습니다.")
        elif err_det:
            print("행렬식 방법으로 역행렬을 구할 수 없습니다.")
        elif err_gj:
            print("가우스-조던 방법으로 역행렬을 구할 수 없습니다.")
        else:
            if compare_matrix(inv_det, inv_gj, size):
                print("두 방법의 결과가 동일합니다.")
                verify_inverse(matrix, inv_det, size)
            else:
                print("두 방법의 결과가 다릅니다.")

    except Exception as e:
        print("예상치 못한 오류가 발생했습니다:", e)


if __name__ == "__main__":
    main()