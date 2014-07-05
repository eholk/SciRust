use matrix::BasicMatrix;

pub fn to_str<T: ToStr, M: BasicMatrix<T>>(m: &M) -> String {
    let mut s = "".to_string();
    for i in range(0, m.num_rows()) {
        s = s + "[ ";
        for j in range(0, m.num_cols()) {
            s = s + m.get(i, j).to_str() + " ";
        }
        s = s + "]\n";
    }
    s
}
