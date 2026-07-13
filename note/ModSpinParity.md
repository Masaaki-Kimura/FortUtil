# ModSpinParity マニュアル

## 概要

**ModSpinParity** は、原子核構造計算などで頻出する **スピンとパリティの組 $(J^\pi)$** を簡潔に扱うためのユーティリティです。
整数表現 $2J$（偶奇判定が容易）を用い、自然パリティや三角条件などの判定も行えます。

本モジュールは主に以下の用途を想定しています：

* `$J^\pi$` を表す構造体として保持（例：`SpinParity(2,POS)` → $J=1^+$）
* 文字列変換（`"2+"`, `"1-"` 形式）
* 三角条件・自然パリティ判定
* クラスター構造・遷移選択則などの補助

より詳細な使用例は **テストコード**（`modspinparity_test.f90`）を参照してください。

---

## 型：`type(SpinParity)`

| 成員                                    | 型                               | 説明 |
| ------------------------------------- | ------------------------------- | -- |
| `integer(int32) :: j`                 | 2J表記のスピン（例：`j=2` → J=1）         |    |
| `integer(int32) :: p`                 | パリティ。`POS`=+, `NEG`=−    |
| `integer(int32), allocatable :: m(:)` | $m=-J, ..., +J$ の配列             |    |
| `integer(int32) :: dim`               | $2J+1$ の次元数                     |    |

### コンストラクタ

#### (1) 文字列から生成

```fortran
type(SpinParity) :: sp
sp = SpinParity("2+")
```

* `"2+"` や `"3-"` のような表記に対応。
* 不正文字列は `error stop`。

#### (2) 整数 + パリティから生成

```fortran
sp = SpinParity(2, POS)  ! J=1, +
sp = SpinParity(3, NEG)  ! J=3/2, -
```

* `j<0` や `p<0` の場合は内部でチェックし、停止。

---

## 型手続き一覧

| 手続き                              | 戻り値            | 概要                         |
| -------------------------------- | -------------- | -------------------------- |
| `function str(this)`             | `character(:)` | `"2+"`, `"1-"` のような文字列を返す  |
| `function dim(this)`             | `integer`      | $2J+1$ を返す                 |
| `function is_natural(this)`      | `logical`      | 自然パリティ判定：$\pi=(-1)^J$      |
| `function is_triangle(this,b,c)` | `logical`      | 三角条件とパリティ保存を同時に判定          |
| `function parity_sign(this)`     | `integer`      | $+1/-1$ の数値を返す             |

---

## 内部動作メモ

* $J$ は内部的に `j = 2J` で保持。
  したがって、`J = j/2`。
* 自然パリティ判定は：

  ```fortran
  is_natural = (mod(j/2, 2) == 0 .and. p == POS) .or. &
                (mod(j/2, 2) == 1 .and. p == NEG)
  ```
* 三角条件は：

  ```fortran
  abs(a.j - b.j) <= c.j .and. (a.j + b.j) >= c.j .and. &
  mod(a.j + b.j + c.j, 2) == 0
  ```

  に加え、パリティ積一致も判定。

---

## 使用例（抜粋）

### 基本操作

```fortran
use ModSpinParity
type(SpinParity) :: a, b, c

a = SpinParity("2+")     ! J=1+
b = SpinParity("2-")     ! J=1−
c = SpinParity(4, POS)   ! J=2+

print *, 'a =', a%str()
print *, 'dim(a) =', a%dim()
```

### 自然パリティの確認

```fortran
a = SpinParity(4, POS)   ! J=2+
if (a%is_natural()) print *, 'Natural parity'
```

### 三角条件の判定

```fortran
a = SpinParity(2, POS)   ! J=1+
b = SpinParity(2, POS)
c = SpinParity(4, POS)
if (a%is_triangle(b, c)) print *, 'Triangle OK'
```

### スピンだけ使いたい場合

```fortran
! ダミーで自然パリティを与えておく
sp = SpinParity(4, POS)  ! J=2+
```

---


## 参考：テスト項目（modspinparity_test.f90）

1. 文字列・整数コンストラクタ
2. `sp%str()` 表示
3. $m$ 配列構築
4. `is_natural()` 判定
5. `is_triangle()` 検証
6. スピンのみ利用例

---

