import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Scanner;

public class replication {

	static double gamma = 1 / 320;
	static double beta = 1;
	static int time = 0;
	static int serviceTimeR = 4;
	static int serviceTimeC = 4;
	static int f = 5;
	static int commitmentOverride = 90;
	static int delta1 = 1;
	static int delta2 = 1;
	static int deltaU = 1;

	// 				0 	1 	2 					3 				4 				5
	// COURIER: 	x 	y 	available (e_d) 	logon (e_c) 	logoff (l_c)
	// ORDER: 		x 	y 	r 					placed (a_o) 	ready (e_o) 	delivered
	// RESTAURANT 	x 	y

	public static void main(String[] args) throws FileNotFoundException {
		File dir = new File("instances");
		File[] folder = dir.listFiles();
		for (File table : folder) {

			File[] filenames = table.listFiles();
			int O = 0;
			int R = 0;
			int C = 0;
			for (File file : filenames) {
				if (file.getName().equals("instance_characteristics.txt")) {
					Scanner characteristics = new Scanner(file);
					String numberOfOrders = characteristics.nextLine();
					O = Integer.valueOf(numberOfOrders.replace("number of orders: ", ""));
					String numberOfRestaurants = characteristics.nextLine();
					R = Integer.valueOf(numberOfRestaurants.replace("number of restaurants: ", ""));
					String numberOfCouriers = characteristics.nextLine();
					C = Integer.valueOf(numberOfCouriers.replace("number of couriers: ", ""));
					characteristics.close();
				}
			}
			int[][] order = new int[O][6];
			int[][] restaurant = new int[R][2];
			int[][] courier = new int[C][5];
			for (File file : filenames) {
				if (file.getName().equals("orders.txt")) {
					Scanner orders = new Scanner(file);

					Scanner orderScanner = new Scanner(orders.nextLine());
					for (int o = 0; o < O; o++) {
						orderScanner = new Scanner(orders.nextLine());
						String index = orderScanner.next();
						order[o][0] = Integer.valueOf(orderScanner.next());
						order[o][1] = Integer.valueOf(orderScanner.next());
						order[o][3] = Integer.valueOf(orderScanner.next());
						order[o][2] = Integer.valueOf(orderScanner.next().replace("r", ""));
						order[o][4] = Integer.valueOf(orderScanner.next());
					}
					orders.close();
				}
				if (file.getName().equals("restaurants.txt")) {
					Scanner restaurants = new Scanner(file);

					Scanner restaurantScanner = new Scanner(restaurants.nextLine());
					for (int r = 0; r < R; r++) {
						restaurantScanner = new Scanner(restaurants.nextLine());
						String index = restaurantScanner.next();
						restaurant[r][0] = Integer.valueOf(restaurantScanner.next());
						restaurant[r][1] = Integer.valueOf(restaurantScanner.next());
					}
					restaurants.close();
				}
				if (file.getName().equals("couriers.txt")) {
					Scanner couriers = new Scanner(file);

					Scanner courierScanner = new Scanner(couriers.nextLine());
					for (int c = 0; c < C; c++) {
						courierScanner = new Scanner(couriers.nextLine());
						String index = courierScanner.next();
						courier[c][0] = Integer.valueOf(courierScanner.next());
						courier[c][1] = Integer.valueOf(courierScanner.next());
						courier[c][2] = Integer.valueOf(courierScanner.next());
						courier[c][3] = courier[c][2];
						courier[c][4] = Integer.valueOf(courierScanner.next());
					}
					couriers.close();
				}
			}
			events(order, restaurant, courier);
		}
	}

	public static void events(int[][] order, int[][] restaurant, int[][] courier) {

		ArrayList<ArrayList<Integer>> S = new ArrayList<>(); // list of active bundles

		while (time < 1440) {
			System.out.println(time);
			ArrayList<Integer> D = new ArrayList<>();

			for (int c = 0; c < courier.length; c++) {
				if (courier[c][3] < time && courier[c][4] > time && courier[c][2] <= time + delta2)
					D.add(c);
			}

			int Z = 0;
			for (int o = 0; o < order.length; o++) {
				if (order[o][4] <= time + deltaU && order[o][5] == 0)
					Z++;
			}

			if (D.size() > 0)
				Z = (int) Math.ceil(Z / D.size());
			Z = Math.max(Z, 1);

			for (int r = 0; r < restaurant.length; r++) {
				// List of upcoming orders
				ArrayList<Integer> U = new ArrayList<>();
				for (int o = 0; o < order.length; o++) {
					if (order[o][2] == r && order[o][4] <= time + deltaU && order[o][5] == 0)
						U.add(o);
				}

				int k = 0;
				for (int c = 0; c < courier.length; c++) {
					if (courier[c][0] == restaurant[r][0] && courier[c][1] == restaurant[r][1])
						k++;
				}
				// System.out.println("Starting procedure 1: " +Z + " " + k + " " + S.size() + "
				// " + U.size());
				S = procedure1(Z, k, S, U, order);

				if (D.size() > 0) {
					// System.out.println("Starting assignment");
					int[][] x = assignment(S, D, r, order, courier, restaurant);
					// System.out.println("Starting commitment");
					x = commitment(S, D, r, x, order, courier, restaurant);
					for (int d : D) {
						for (ArrayList<Integer> s : S) {
							if (x[S.indexOf(s)][D.indexOf(d)] == 1) {
								courier[d][0] = order[s.get(s.size() - 1)][0];
								courier[d][1] = order[s.get(s.size() - 1)][1];
								courier[d][2] += travelTime(d, r, courier, restaurant) + serviceTimeR
										+ travelTime(d, s.get(0), courier, order) + travelTime(s, order);
								for (int o : s) {
									order[o][5] = 1;
								}
							}
							if (x[S.indexOf(s)][D.indexOf(d)] == 1 / 2) {
								courier[d][0] = order[s.get(s.size() - 1)][0];
								courier[d][1] = order[s.get(s.size() - 1)][1];
								courier[d][2] += travelTime(d, r, courier, restaurant);
							}
						}
					}
				}

			}
			time += f;
		}
	}

	public static ArrayList<ArrayList<Integer>> procedure1(int Z, int k, ArrayList<ArrayList<Integer>> S,
			ArrayList<Integer> U, int[][] order) {

		U = sortOrders(U, order); // sorted set of upcoming orders

		int m = (int) Math.max(k, Math.ceil(U.size() / Z)); // target number of bundles

		for (int i = S.size(); i < m; i++) {
			S.add(new ArrayList<>());
		}

		for (int u : U) {
			boolean ok = false;

			int s = bestInsert(S, order, u)[0];
			int i = bestInsert(S, order, u)[1];

			boolean e = (bestInsert(S, order, u)[2] != 0);
			if (S.get(s).size() < Z || e) {
				S.get(s).add(i, u);

			}
		}
		// ADD REMOVE-REINSERT
		return S;
	}

	public static ArrayList<Integer> sortOrders(ArrayList<Integer> U, int[][] order) {
		for (int u = 1; u < U.size(); u++) {
			int o = U.get(u);
			if (order[o][4] < order[U.get(u - 1)][4]) {
				U.remove(U.indexOf(o));
				for (int i = 0; i <= u; i++) {
					if (order[o][4] < order[U.get(i)][4]) {
						U.add(i, o);
						break;
					}
				}
			}
		}
		return U;
	}

	/**
	 * 
	 * @param S     set of bundles from restaurant r
	 * @param order set of orders
	 * @param u     order to be inserted
	 * @return best route, best position, efficiency
	 */
	public static int[] bestInsert(ArrayList<ArrayList<Integer>> S, int[][] order, int u) {
		int route = 0;
		int pos = 0;

		double cost = 0;
		double minCost = Integer.MAX_VALUE;

		for (ArrayList<Integer> s : S) {

			for (int o : s) {
				cost = travelTime(u, s, s.indexOf(o), order) + delayCost(u, s, s.indexOf(o), order);
				if (cost < minCost) {
					minCost = cost;
					route = S.indexOf(s);
					pos = s.indexOf(o);
				}
			}
		}

		int efficiency = (travelTime(S.get(route), order)
				/ S.get(route).size() > travelTime(u, S.get(route), pos, order) / (S.get(route).size() + 1)) ? 1 : 0;

		// = 1 if time per order decreases with insertion

		return new int[] { route, pos, efficiency };
	}

	public static double delayCost(int u, ArrayList<Integer> s, int o, int[][] order) {
		ArrayList<Integer> copyS = new ArrayList<>(s);
		copyS.add(o, u);
		return delayCost(copyS, order);
	}

	public static double delayCost(ArrayList<Integer> s, int[][] order) {
		double C = 0; // total travel time
		int maxE = 0;

		for (int o : s) {
			maxE = Math.max(maxE, order[o][4]);
		}

		for (int o : s) {
			C += beta * (maxE - order[o][4]);
		}
		return C;
	}

	/**
	 * Travel time/cost of inserting order u into bundle s at position o
	 * 
	 * @param u     the order
	 * @param s     the bundle
	 * @param o     the position
	 * @param order set of orders
	 * @return Euclidean distance with gamma speed multiplier
	 */

	public static double travelTime(int u, ArrayList<Integer> s, int o, int[][] order) {
		ArrayList<Integer> copyS = new ArrayList<>(s);
		copyS.add(o, u);
		return travelTime(copyS, order);
	}

	/**
	 * Total travel time/cost through all customers in a bundle;
	 * 
	 * @param s     the bundle
	 * @param order set of orders
	 * @param b     service delay cost
	 * @return Euclidean distance with gamma speed multiplier
	 */
	public static double travelTime(ArrayList<Integer> s, int[][] order) {

		double T = serviceTimeC; // total travel time

		for (int i = 0; i < s.size() - 1; i++) {
			T += travelTime(i, i + 1, order, order);
			T += serviceTimeC; // service time at i+1;
		}
		return T;
	}

	/**
	 * Travel time from point A of type a to point B of type b
	 * 
	 * @param A
	 * @param B
	 * @param a
	 * @param b
	 * @return the travel time, with speed gamma
	 */
	public static double travelTime(int A, int B, int[][] a, int[][] b) {

		int x1 = a[A][0];
		int x2 = b[B][0];
		int y1 = a[A][1];
		int y2 = b[B][1];

		return gamma * Math.ceil(Math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)));
	}

	public static int[][] assignment(ArrayList<ArrayList<Integer>> S, ArrayList<Integer> D, int r, int[][] order,
			int[][] courier, int[][] restaurant) {
		int[][] x = new int[S.size()][D.size()];
		int[] assigned = new int[D.size()];
		for (ArrayList<Integer> s : S) {
			double objective = 0;
			double bestObjective = Integer.MIN_VALUE;
			int bestCourier = D.get(0);
			int N = s.size();
			double theta = beta;
			for (int d : D) {
				int maxE = Integer.MIN_VALUE;
				for (int o : s) {
					maxE = Math.max(maxE, order[o][5]);
				}
				double pi = Math.max(maxE, time + travelTime(r, d, restaurant, courier) + 1 / 2 * serviceTimeR); // bundle
				double maxDelta = pi + 1 / 2 * serviceTimeR + travelTime(r, s.get(0), restaurant, order)
						+ travelTime(s, order); // dropoff time // pickup
				objective = (N / maxDelta - courier[d][3]) - theta * (pi - maxE);
				if (objective > bestObjective && assigned[D.indexOf(d)] == 0) {
					bestObjective = objective;
					bestCourier = d;
				}
			}

			x[S.indexOf(s)][D.indexOf(bestCourier)] = 1;
			assigned[D.indexOf(bestCourier)] = 1;
		}

		return x;
	}

	public static int[][] commitment(ArrayList<ArrayList<Integer>> S, ArrayList<Integer> D, int r, int[][] x,
			int[][] order, int[][] courier, int[][] restaurant) {
		for (ArrayList<Integer> s : S) {
			for (int d : D) {
				boolean reach = courier[d][3] + travelTime(d, r, courier, restaurant) < time + f;
				int maxE = Integer.MIN_VALUE;
				for (int o : s) {
					maxE = Math.max(maxE, order[o][5]);
				}
				boolean ready = maxE < time + f;
				boolean start = courier[d][3] < time + f;
				int minA = Integer.MAX_VALUE;
				for (int o : s) {
					maxE = Math.min(minA, order[o][4]);
				}
				boolean late = minA + commitmentOverride < time;
				if (!late && !start) {
					x[S.indexOf(s)][D.indexOf(d)] = 0;
				} else if (!late && (!reach || !ready)) {
					x[S.indexOf(s)][D.indexOf(d)] = 1 / 2;
				}
			}
		}
		return x;
	}

}
