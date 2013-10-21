use std::local_data;

use matrix::BasicMatrix;

pub fn to_str<T: ToStr, M: BasicMatrix<T>>(m: &M) -> ~str {
    let mut s = ~"";
    for i in range(0, m.num_rows()) {
        s = s + "[ ";
        for j in range(0, m.num_cols()) {
            s = s + m.get(i, j).to_str() + " ";
        }
        s = s + "]\n";
    }
    s
}

local_data_key!(CUR_INDENT: uint)

fn indent() -> ~str {
    "  ".repeat(local_data::get(CUR_INDENT,
                                |opt| match opt {
                                    Some(x) => *x,
                                    None => 0
                                }))
}

pub fn tracefn<T, U: ToStr>(m: U, f: &fn() -> T) -> T {
    return f();

    let m = m.to_str();
    info!("{:s}", indent() + m);
    local_data::modify(CUR_INDENT,
                       |opt| match opt {
                           Some(x) => Some(x + 1),
                           None => Some(1)
                       });
    
    let r = f();

    local_data::modify(CUR_INDENT,
                       |opt| match opt {
                           Some(x) => Some(x - 1),
                           None => fail!("indent level too low")
                       });
    info!("{:s}", (indent() + "~" + m));
    
    r
}
