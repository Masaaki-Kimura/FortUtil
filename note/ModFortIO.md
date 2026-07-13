# ModFortIO マニュアル

## 概要
**ModFortIO** は Fortran での I/O を安全・実用的に行うためのユーティリティ集です。  
- 端末の**カラー出力**  
- **ファイル名操作**（パス・拡張子・ベース名など）  
- **ファイル存在確認**  
- **NaN 生成**  
- **コマンドライン引数処理**（位置指定・逐次取得）  
- **`FileHandler` クラス**による安全な open/close ラッパ  

より詳細な使い方は、付属の **テストコード**（`test_modfortio.f90`）を参照してください。

---

## 公開定数・変数
- `integer(int32), parameter :: MAXBUF = 512`  
  汎用バッファ長
- `integer(int32), parameter :: STDIN, STDOUT, STDERR`  
  `use iso_fortran_env, only: input_unit, output_unit, error_unit` に基づく予約ユニット
- `character(len=1), parameter :: ASCII_TAB, ASCII_LF, ASCII_CR, ASCII_ESC`
- `logical, save :: modfortio_enable_text_color = .true.`  
  端末色の有効/無効
- `character(len=1) :: modfortio_error_color = 'r'`  
  既定のエラーカラー（`'r','g','y','b','m','c'`）

---

## 型：`type(FileHandler)`
| 成員 | 説明 |
|---|---|
| `integer(int32) :: unit` | 接続ユニット番号。初期値は `STDERR` |
| `character(:), allocatable :: name` | 接続中のファイル名 |

### 型手続き
- `logical function open(this, fn, status, action, position, form, msg)`  
  - 既定：`status='unknown'`, `action='readwrite'`, `position='rewind'`, `form='formatted'`  
  - `form='unformatted'` の場合は `access='stream'` を強制  
  - 既に同名ファイルが開かれている場合はエラー（`false` を返す）  
  - 失敗時 `msg`（可変長文字列）に詳細
- `logical function open_binary(this, fn, status, action, position, msg)`  
  - `open(..., form='unformatted')` の薄いラッパ
- `logical function close(this)`  
  - 正常時 `true`、未接続なら警告出力して `false`  
  - 成功時に `unit=STDERR` へ戻す
- `logical function create_blank_file(this, fn)`  
  - `status='new'` で空ファイルを作成してすぐ close  
  - 既存なら `false`


---

## 関数・サブルーチン

### 端末カラー
- `pure function text_color(color, src) result(s)`  
  - `modfortio_enable_text_color=.false.` のときは `src` をそのまま返す  
  - `color∈{'r','g','y','b','m','c'}`

### パス操作
- `pure function append_slash(src)`
- `pure function file_path(src)`  
  - 例：`'/a/b/c.txt' -> '/a/b/'`、スラッシュ無しなら `'./'`
- `pure function file_name(src)`  
  - 例：`'/a/b/c.txt' -> 'c.txt'`
- `pure function file_basename(src)`  
  - 例：`'/a/b/c.tar.gz' -> 'c.tar'`
- `pure function file_ext(src)`  
  - 例：`'/a/b/c.tar.gz' -> '.gz'`、拡張子がなければ空文字

### ファイル存在確認
- `logical function file_exists(fn)`

### 数値ユーティリティ
- `pure real(real64) function ieee_nan()`  
  - `ieee_value(0.0_real64, ieee_quiet_nan)`

### コマンドライン
- `function get_arg(n) result(s)`  
  - 第 `n` 引数を返す（なければ空文字）
- `function next_arg() result(s)`  
  - 呼ぶたびに次引数を返す（内部 `save` カーソル）  
  - **注**：本関数は引数を取りません

---

## 使用例（抜粋）

### 基本的なファイル書き出し

```fortran
use ModFortIO
type(FileHandler) :: fh
logical :: ok
character(:), allocatable :: msg

! テキストファイルを新規作成（上書き可）
ok = fh%open('out.txt', status='replace', action='write', form='formatted', msg=msg)
if (.not. ok) error stop trim(msg)

! 出力例
write(fh%unit,'(A)') 'hello world'
write(fh%unit,'(A)') 'this is ModFortIO test.'

! クローズ
ok = fh%close()
if (.not. ok) error stop 'close failed'
```

### ファイル存在チェックと空ファイル作成
```fortran
if (.not. file_exists('newfile.dat')) then
    call fh%create_blank_file('newfile.dat')
    print *, 'Created: newfile.dat'
else
    print *, 'Already exists.'
end if
```

### STREAM バイナリ書き込み
```fortran
integer(int32) :: value = 42
ok = fh%open_binary('data.bin', status='replace', action='write', msg=msg)
if (.not. ok) error stop trim(msg)
write(fh%unit) value
ok = fh%close()
```

### パス・ファイル名操作
```fortran
print *, file_path('/usr/local/test.txt')     ! => '/usr/local/'
print *, file_name('/usr/local/test.txt')     ! => 'test.txt'
print *, file_basename('/usr/local/test.txt') ! => 'test'
print *, file_ext('/usr/local/test.txt')      ! => '.txt'
```

### コマンドライン引数の取得
```fortran
character(:), allocatable :: arg
arg = get_arg(1)
print *, 'arg1 = ', trim(arg)

arg = next_arg()
print *, 'next = ', trim(arg)
```

### NaN 生成と確認
```fortran
real(real64) :: x
x = ieee_nan()
if (x == x) error stop 'ieee_nan failed'  ! NaNは自身と等しくない
```

### カラー出力（端末依存）
```fortran
use ModFortIO
print *, text_color('g', 'Success: file created')
print *, text_color('r', 'Error: invalid argument')
```
