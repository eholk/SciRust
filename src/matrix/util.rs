use matrix::BasicMatrix;

pub fn to_str<T: Copy ToStr, M: BasicMatrix<T>>(m: &M) -> ~str {
    let mut s = ~"";
    for uint::range(0, m.num_rows()) |i| {
        s += ~"[ ";
        for uint::range(0, m.num_cols()) |j| {
            s += m.get(i, j).to_str() + " ";
        }
        s += "]\n";
    }
    move s
}
