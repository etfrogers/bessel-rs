use std::error::Error;

use crate::{bessel_j, bessel_y};

use crate::{BesselFunType, bessel_zeros};
use approx::assert_relative_eq;
use rstest::rstest;

const J_ZEROS: &str = "#
#k	J_0(x)	J_1(x)	J_2(x)	J_3(x)	J_4(x)	J_5(x)
1	2.4048	3.8317	5.1356	6.3802	7.5883	8.7715
2	5.5201	7.0156	8.4172	9.7610	11.0647	12.3386
3	8.6537	10.1735	11.6198	13.0152	14.3725	15.7002
4	11.7915	13.3237	14.7960	16.2235	17.6160	18.9801
5	14.9309	16.4706	17.9598	19.4094	20.8269	22.2178";

const JP_ZEROS: &str = "#
#k	J_0'(x)	J_1'(x)	J_2'(x)	J_3'(x)	J_4'(x)	J_5'(x)
1	0.0000	1.8412	3.0542	4.2012	5.3175	6.4156
2	3.8317	5.3314	6.7061	8.0152	9.2824	10.5199
3	7.0156	8.5363	9.9695	11.3459	12.6819	13.9872
4	10.1735	11.7060	13.1704	14.5858	15.9641	17.3128
5	13.3237	14.8636	16.3475	17.7887	19.1960	20.5755";

//	16.4706

const Y_ZEROS: &str = "#
0.89357697	3.95767842	7.08605106	10.22234504	13.36109747
2.19714133	5.42968104	8.59600587	11.74915483	14.89744213
3.38424177	6.79380751	10.02347798	13.20998671	16.37896656
4.52702466	8.09755376	11.39646674	14.62307774	17.81845523
5.64514789	9.36162062	12.73014447	15.99962709	19.22442896";

const YP_ZEROS: &str = "#
2.19714133	5.42968104	8.59600587	11.74915483	14.89744213
3.68302286	6.94149995	10.12340466	13.28575816	16.44005801
5.00258293	8.3507247	11.57419547	14.76090931	17.93128594
6.25363321	9.69878798	12.97240905	16.1904472	19.38238845
7.46492174	11.00516915	14.33172352	17.58443602	20.80106234";

fn get_records(input: &str) -> Result<Vec<Vec<String>>, Box<dyn Error>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .has_headers(false)
        .from_reader(input.as_bytes());
    Ok(rdr
        .records()
        .map(|r| r.unwrap().iter().map(|s| s.to_string()).collect())
        .collect())
}

fn parse_zeros_python(input: &str) -> Vec<Vec<f64>> {
    get_records(input)
        .unwrap()
        .into_iter()
        .map(|line| line.iter().map(|v| v.parse().unwrap()).collect())
        .collect()
}

fn parse_zeros_wa(input: &str) -> Vec<Vec<f64>> {
    let vv: Vec<Vec<f64>> = get_records(input)
        .unwrap()
        .into_iter()
        .map(|line| line.iter().map(|v| v.parse().unwrap()).collect())
        .collect();
    let mut arr = vec![vec![]; vv[0].len()];
    for row in &vv {
        for (i, val) in row[1..].iter().enumerate() {
            arr[i].push(*val);
        }
    }
    arr
}

#[rstest]
#[case(BesselFunType::J, parse_zeros_wa(J_ZEROS))]
#[case(BesselFunType::JP, parse_zeros_wa(JP_ZEROS))]
#[case(BesselFunType::Y, parse_zeros_python(Y_ZEROS))]
#[case(BesselFunType::YP, parse_zeros_python(YP_ZEROS))]
fn test_against_hard_coded_zeros(#[case] fun_type: BesselFunType, #[case] zeros: Vec<Vec<f64>>) {
    for (order, expected) in (0_i32..).zip(zeros) {
        let actual = bessel_zeros(&fun_type, order, expected.len(), 1e-6);
        // println!("Order {}", order);
        // println!("{:?}", actual);
        // println!("{:?}", expected);
        expected
            .into_iter()
            .zip(actual)
            .for_each(|(ev, av)| assert_relative_eq!(ev, av, epsilon = 5e-4));
    }
}

#[rstest]
fn test_evaluation_at_zero_j() {
    for order in 0..20 {
        println!("{}", order);
        let zeros = bessel_zeros(&BesselFunType::J, order, 100, 0.1e-6);
        for v in zeros {
            print!("{}, ", v);
            assert_relative_eq!(0.0, bessel_j(order, v).unwrap(), epsilon = 1e-6)
        }
    }
}

#[rstest]
fn test_evaluation_at_zero_y() {
    for order in 0..20 {
        let zeros = bessel_zeros(&BesselFunType::Y, order, 100, 0.1e-6);
        for v in zeros {
            assert_relative_eq!(0.0, bessel_y(order, v).unwrap(), epsilon = 1e-6);
        }
    }
}
