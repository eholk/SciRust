use matrix::BasicMatrix;

pub fn to_str<T: ToStr, M: BasicMatrix<T>>(m: &M) -> ~str {
    let mut s = "".to_owned();
    for i in range(0, m.num_rows()) {
        s = s + "[ ";
        for j in range(0, m.num_cols()) {
            s = s + m.get(i, j).to_str() + " ";
        }
        s = s + "]\n";
    }
    s
}
