using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Diagnostics;


namespace TSP
{

    class ProblemAndSolver
    {

        private class TSPSolution
        {
            /// <summary>
            /// we use the representation [cityB,cityA,cityC] 
            /// to mean that cityB is the first city in the solution, cityA is the second, cityC is the third 
            /// and the edge from cityC to cityB is the final edge in the path.  
            /// You are, of course, free to use a different representation if it would be more convenient or efficient 
            /// for your data structure(s) and search algorithm. 
            /// </summary>
            public ArrayList
                Route;

            /// <summary>
            /// constructor
            /// </summary>
            /// <param name="iroute">a (hopefully) valid tour</param>
            public TSPSolution(ArrayList iroute)
            {
                Route = new ArrayList(iroute);
            }

            /// <summary>
            /// Compute the cost of the current route.  
            /// Note: This does not check that the route is complete.
            /// It assumes that the route passes from the last city back to the first city. 
            /// </summary>
            /// <returns></returns>
            public double costOfRoute()
            {
                // go through each edge in the route and add up the cost. 
                int x;
                City here;
                double cost = 0D;

                for (x = 0; x < Route.Count - 1; x++)
                {
                    here = Route[x] as City;
                    cost += here.costToGetTo(Route[x + 1] as City);
                }

                return cost;
            }
        }

        #region Private members 

        /// <summary>
        /// Default number of cities (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Problem Size text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int DEFAULT_SIZE = 25;

        /// <summary>
        /// Default time limit (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Time text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int TIME_LIMIT = 60;        //in seconds

        private const int CITY_ICON_SIZE = 5;


        // For normal and hard modes:
        // hard mode only
        private const double FRACTION_OF_PATHS_TO_REMOVE = 0.20;

        /// <summary>
        /// the cities in the current problem.
        /// </summary>
        private City[] Cities;
        /// <summary>
        /// a route through the current problem, useful as a temporary variable. 
        /// </summary>
        private ArrayList Route;
        /// <summary>
        /// best solution so far. 
        /// </summary>
        private TSPSolution bssf; 

        /// <summary>
        /// how to color various things. 
        /// </summary>
        private Brush cityBrushStartStyle;
        private Brush cityBrushStyle;
        private Pen routePenStyle;


        /// <summary>
        /// keep track of the seed value so that the same sequence of problems can be 
        /// regenerated next time the generator is run. 
        /// </summary>
        private int _seed;
        /// <summary>
        /// number of cities to include in a problem. 
        /// </summary>
        private int _size;

        /// <summary>
        /// Difficulty level
        /// </summary>
        private HardMode.Modes _mode;

        /// <summary>
        /// random number generator. 
        /// </summary>
        private Random rnd;

        /// <summary>
        /// time limit in milliseconds for state space search
        /// can be used by any solver method to truncate the search and return the BSSF
        /// </summary>
        private int time_limit;
        #endregion

        #region Public members

        /// <summary>
        /// These three constants are used for convenience/clarity in populating and accessing the results array that is passed back to the calling Form
        /// </summary>
        public const int COST = 0;           
        public const int TIME = 1;
        public const int COUNT = 2;
        
        public int Size
        {
            get { return _size; }
        }

        public int Seed
        {
            get { return _seed; }
        }
        #endregion

        #region Constructors
        public ProblemAndSolver()
        {
            this._seed = 1; 
            rnd = new Random(1);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed)
        {
            this._seed = seed;
            rnd = new Random(seed);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed, int size)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = TIME_LIMIT * 1000;                        // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        public ProblemAndSolver(int seed, int size, int time)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = time*1000;                        // time is entered in the GUI in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        #endregion

        #region Private Methods

        /// <summary>
        /// Reset the problem instance.
        /// </summary>
        private void resetData()
        {

            Cities = new City[_size];
            Route = new ArrayList(_size);
            bssf = null;

            if (_mode == HardMode.Modes.Easy)
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble());
            }
            else // Medium and hard
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble(), rnd.NextDouble() * City.MAX_ELEVATION);
            }

            HardMode mm = new HardMode(this._mode, this.rnd, Cities);
            if (_mode == HardMode.Modes.Hard)
            {
                int edgesToRemove = (int)(_size * FRACTION_OF_PATHS_TO_REMOVE);
                mm.removePaths(edgesToRemove);
            }
            City.setModeManager(mm);

            cityBrushStyle = new SolidBrush(Color.Black);
            cityBrushStartStyle = new SolidBrush(Color.Red);
            routePenStyle = new Pen(Color.Blue,1);
            routePenStyle.DashStyle = System.Drawing.Drawing2D.DashStyle.Solid;
        }

        #endregion

        #region Public Methods

        /// <summary>
        /// make a new problem with the given size.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode)
        {
            this._size = size;
            this._mode = mode;
            resetData();
        }

        /// <summary>
        /// make a new problem with the given size, now including timelimit paremeter that was added to form.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode, int timelimit)
        {
            this._size = size;
            this._mode = mode;
            this.time_limit = timelimit*1000;                                   //convert seconds to milliseconds
            resetData();
        }

        /// <summary>
        /// return a copy of the cities in this problem. 
        /// </summary>
        /// <returns>array of cities</returns>
        public City[] GetCities()
        {
            City[] retCities = new City[Cities.Length];
            Array.Copy(Cities, retCities, Cities.Length);
            return retCities;
        }

        /// <summary>
        /// draw the cities in the problem.  if the bssf member is defined, then
        /// draw that too. 
        /// </summary>
        /// <param name="g">where to draw the stuff</param>
        public void Draw(Graphics g)
        {
            float width  = g.VisibleClipBounds.Width-45F;
            float height = g.VisibleClipBounds.Height-45F;
            Font labelFont = new Font("Arial", 10);

            // Draw lines
            if (bssf != null)
            {
                // make a list of points. 
                Point[] ps = new Point[bssf.Route.Count];
                int index = 0;
                foreach (City c in bssf.Route)
                {
                    if (index < bssf.Route.Count -1)
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[index+1]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    else 
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[0]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    ps[index++] = new Point((int)(c.X * width) + CITY_ICON_SIZE / 2, (int)(c.Y * height) + CITY_ICON_SIZE / 2);
                }

                if (ps.Length > 0)
                {
                    g.DrawLines(routePenStyle, ps);
                    g.FillEllipse(cityBrushStartStyle, (float)Cities[0].X * width - 1, (float)Cities[0].Y * height - 1, CITY_ICON_SIZE + 2, CITY_ICON_SIZE + 2);
                }

                // draw the last line. 
                g.DrawLine(routePenStyle, ps[0], ps[ps.Length - 1]);
            }

            // Draw city dots
            foreach (City c in Cities)
            {
                g.FillEllipse(cityBrushStyle, (float)c.X * width, (float)c.Y * height, CITY_ICON_SIZE, CITY_ICON_SIZE);
            }

        }

        /// <summary>
        ///  return the cost of the best solution so far. 
        /// </summary>
        /// <returns></returns>
        public double costOfBssf ()
        {
            if (bssf != null)
                return (bssf.costOfRoute());
            else
                return -1D; 
        }

        /// <summary>
        /// This is the entry point for the default solver
        /// which just finds a valid random tour 
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] defaultSolveProblem()
        {
            int i, swap, temp, count=0;
            string[] results = new string[3];
            int[] perm = new int[Cities.Length];
            Route = new ArrayList();
            Random rnd = new Random();
            Stopwatch timer = new Stopwatch();

            timer.Start();

            do
            {
                for (i = 0; i < perm.Length; i++)                                 // create a random permutation template
                    perm[i] = i;
                for (i = 0; i < perm.Length; i++)
                {
                    swap = i;
                    while (swap == i)
                        swap = rnd.Next(0, Cities.Length);
                    temp = perm[i];
                    perm[i] = perm[swap];
                    perm[swap] = temp;
                }
                Route.Clear();
                for (i = 0; i < Cities.Length; i++)                            // Now build the route using the random permutation 
                {
                    Route.Add(Cities[perm[i]]);
                }
                bssf = new TSPSolution(Route);
                count++;
            } while (costOfBssf() == double.PositiveInfinity);                // until a valid route is found
            timer.Stop();

            results[COST] = costOfBssf().ToString();                          // load results array
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString();

            return results;
        }

        public ArrayList getRandomRoute()
        {
            int i, swap, temp = 0;
            int[] perm = new int[Cities.Length];
            Route = new ArrayList();
            Random rnd = new Random();


            for (i = 0; i < perm.Length; i++)                                 // create a random permutation template
                perm[i] = i;
            for (i = 0; i < perm.Length; i++)
            {
                swap = i;
                while (swap == i)
                    swap = rnd.Next(0, Cities.Length);
                temp = perm[i];
                perm[i] = perm[swap];
                perm[swap] = temp;
            }
            Route.Clear();
            for (i = 0; i < Cities.Length; i++)                            // Now build the route using the random permutation 
            {
                Route.Add(Cities[perm[i]]);
            }

            return Route;
        }

        /// <summary>
        /// performs a Branch and Bound search of the state space of partial tours
        /// stops when time limit expires and uses BSSF as solution
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] bBSolveProblem()
        {
            //Start a timer
            Stopwatch timer = new Stopwatch();
            timer.Start();

            //Initialize variables
            string[] results = new string[3];
            bool found_optimal = false;
            int count = 0;
            int max_queue_size = 0;
            int state_count = 0;
            int bssf_update_count = 0;
            int pruned_count = 0;

            //Get initial greedy path
            greedySolveProblem();
            ArrayList route = bssf.Route;

            double avg_path = bssf.costOfRoute() / Cities.Length;

            //Create initial state and add it to priority queue
            State s = new State(avg_path);
            s.cost_matrix = new double[Cities.Length, Cities.Length];
            for (int row = 0; row < Cities.Length; row++)
            {
                for (int col = 0; col < Cities.Length; col++)
                {
                    s.cost_matrix[row, col] = Cities[row].costToGetTo(Cities[col]);

                    if (row == col)
                    {
                        s.cost_matrix[row, col] = double.PositiveInfinity;
                    }
                }
            }
            s.path.Add(0);
            s.reduceMatrix();
            PriorityQueue<State> q = new PriorityQueue<State>();
            q.add(s);

            //Loop until out of time or the optimal solution is found
            while (timer.ElapsedMilliseconds < time_limit && !found_optimal)
            {
                max_queue_size = Math.Max(max_queue_size, q.Count());

                s = q.remove();
                state_count++;
                if (s.lower_bound < costOfBssf())
                {
                    //Look at children

                    if(s.path.Count == Cities.Length+1)
                    {
                        route = new ArrayList();
                        for(int i=0; i<s.path.Count; i++)
                        {
                            route.Add((City) Cities[(int) s.path[i]]);
                        }
                        TSPSolution tsp = new TSPSolution(route);
                        if(tsp.costOfRoute() < costOfBssf())
                        {
                            bssf = tsp;
                            bssf_update_count++;
                        }
                        else
                        {
                            pruned_count++;
                        }
                        count++;
                    }

                    else if (s.path.Count == Cities.Length)
                    {
                        //Last one, check if you can get back to the beginning
                        q.add(s.getNextState(0));
                    }
                    else
                    {
                        //Not ready for a full path yet

                        for (int i = 0; i < Cities.Length; i++)
                        {
                            if (!s.path.Contains(i))
                            {
                                //Explore
                                q.add(s.getNextState(i));
                                //new_state.reduceMatrix();
                            }
                        }

                    }
                }
                else
                {
                    pruned_count++;
                }

                if (q.Count() == 0)
                {
                    found_optimal = true;
                }

            }

            //Stop the timer
            timer.Stop();

            //Set return values
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString();
            results[COST] = costOfBssf().ToString();

            //Output table data to console
            string table_data = Cities.Length.ToString() + " " + Seed.ToString() + " " + timer.ElapsedMilliseconds/1000.0 + " " + results[COST];
            if (found_optimal)
            {
                table_data = table_data + "*";
            }

            table_data = table_data + " " + max_queue_size.ToString() + " " + bssf_update_count + " " + state_count + " " + pruned_count;
            Debug.Print(table_data);

            return results;
        }

        /////////////////////////////////////////////////////////////////////////////////////////////
        // These additional solver methods will be implemented as part of the group project.
        ////////////////////////////////////////////////////////////////////////////////////////////

        /// <summary>
        /// finds the greedy tour starting from each city and keeps the best (valid) one
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] greedySolveProblem()
        {
            Stopwatch timer = new Stopwatch();
            timer.Start();

            string[] results = new string[3];

            ArrayList route = new ArrayList(); 
            List<int> visited_cities = new List<int>();
            for (int i = 0; i < Cities.Length; i++)
            {
                if (i == 0)
                {
                    route.Add(Cities[i]);
                    visited_cities.Add(i);
                }
                else
                {
                    double shortest_distance = double.PositiveInfinity;
                    int best_j = -1;
                    for (int j = 0; j < Cities.Length; j++)
                    {
                        if (!visited_cities.Contains(j))
                        {
                            City c = route[i - 1] as City;
                            double new_distance = c.costToGetTo(Cities[j]);
                            if (new_distance < shortest_distance)
                            {
                                shortest_distance = new_distance;
                                best_j = j;
                            }
                        }
                    }
                    route.Add(Cities[best_j]);
                    visited_cities.Add(best_j);
                }

            }
            TSPSolution new_solution = new TSPSolution(route);
            bssf = new_solution;

            //TODO: Check if the solution is valid

            timer.Stop();

            results[COST] = costOfBssf().ToString();    // load results into array here, replacing these dummy values
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = "1";

            return results;
        }

        public string[] fancySolveProblem()
        {
            string[] results = new string[3];

            // TODO: Add your implementation for your advanced solver here.

            results[COST] = "not implemented";    // load results into array here, replacing these dummy values
            results[TIME] = "-1";
            results[COUNT] = "-1";

            return results;
        }
        #endregion
    }

    public class State : IComparable<State>
    {
        public ArrayList path = new ArrayList();
        public double lower_bound = 0;
        public double[,] cost_matrix;
        double avg_path;

        public State(double avg_path)
        {
            this.avg_path = avg_path;
        }

        public int CompareTo(State other)
        {
            double self_value = lower_bound - (avg_path * Math.Exp(path.Count));
            double other_value = other.lower_bound - (avg_path * Math.Exp(other.path.Count));

            if (self_value < other_value)
            {
                return -1;
            }
            else if (self_value > other_value)
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }

        public void printCostMatrix()
        {
            string output = "";
            for(var i=0; i<cost_matrix.GetLength(0); i++)
            {
                for(var j=0; j<cost_matrix.GetLength(1); j++)
                {
                    output = output + cost_matrix[i, j].ToString() + " ";
                }
                output = output + "\n";
            }

            Debug.Print(output);
        }

        public void reduceMatrix()
        {
            for(int i=0; i<cost_matrix.GetLength(0); i++)
            {
                double additional_cost = reduceRow(i);
                if (additional_cost != double.PositiveInfinity)
                {
                    lower_bound = lower_bound + additional_cost;
                }
            }

            for(int i=0; i<cost_matrix.GetLength(1); i++)
            {
                double additional_cost = reduceCol(i);
                if (additional_cost != double.PositiveInfinity)
                {
                    lower_bound = lower_bound + additional_cost;
                }
            }
        }

        private double reduceRow(int row_id)
        {
            double lowest_value = double.PositiveInfinity;
            for(int i=0; i<cost_matrix.GetLength(1); i++)
            {
                if (cost_matrix[row_id, i] < lowest_value)
                {
                    lowest_value = cost_matrix[row_id, i];
                }
            }
            for(int i=0; i<cost_matrix.GetLength(1); i++)
            {
                cost_matrix[row_id, i] = cost_matrix[row_id, i] - lowest_value;
            }

            return lowest_value;
        }

        private double reduceCol(int col_id)
        {
            double lowest_value = double.PositiveInfinity;
            for (int i = 0; i < cost_matrix.GetLength(0); i++)
            {
                if (cost_matrix[i, col_id] < lowest_value)
                {
                    lowest_value = cost_matrix[i, col_id];
                }
            }
            for (int i = 0; i < cost_matrix.GetLength(0); i++)
            {
                cost_matrix[i,col_id] = cost_matrix[i,col_id] - lowest_value;
            }

            return lowest_value;
        }

        public State getNextState(int to_idx)
        {
            State s = new State(avg_path);
            s.lower_bound = lower_bound;
            s.path = (ArrayList) path.Clone();
            s.cost_matrix = (double[,])cost_matrix.Clone();

            int from_idx = (int) path[path.Count -1];

            s.lower_bound = s.lower_bound + s.cost_matrix[from_idx, to_idx];

            //Wipe out stuff
            s.cost_matrix[to_idx, from_idx] = double.PositiveInfinity;
            for(int i=0; i<s.cost_matrix.GetLength(0); i++)
            {
                s.cost_matrix[from_idx, i] = double.PositiveInfinity;
                s.cost_matrix[i, to_idx] = double.PositiveInfinity;
            }


            s.path.Add(to_idx);
            s.reduceMatrix();

            return s;

        }
    }

    public class PriorityQueue<T> where T : IComparable<T>
    {
        private List<T> data;

        public PriorityQueue()
        {
            this.data = new List<T>();
        }

        public void add(T item)
        {
            data.Add(item);

            //Initialize the child index to the end of the list
            int child_idx = data.Count - 1;

            while (child_idx > 0)
            {
                int parent_idx = (child_idx - 1) / 2;

                //If the child is smaller than its parent 
                if (data[child_idx].CompareTo(data[parent_idx]) < 0)
                {
                    //Switch the child and parent
                    T tmp = data[child_idx];
                    data[child_idx] = data[parent_idx];
                    data[parent_idx] = tmp;
                    child_idx = parent_idx;
                }
                else
                {
                    //The child is larger than or equal to its parent
                    break;
                }
            }
        }

        public T remove()
        {
            int last_idx = data.Count - 1;
            T front_item = data[0];

            //Move the bottom item to the top
            data[0] = data[last_idx];
            data.RemoveAt(last_idx);
            last_idx--;

            int parent_idx = 0;
            while (true)
            {
                int left_child = parent_idx * 2 + 1;
                
                //If there are no children break
                if(left_child > last_idx)
                {
                    break;
                }

                int right_child = left_child + 1;

                //If the right child is smaller than the left child, use the right child
                if(right_child <= last_idx && data[right_child].CompareTo(data[left_child]) < 0)
                {
                    left_child = right_child;
                }

                //Break if the parent is smaller or equal to the smallest child
                if (data[parent_idx].CompareTo(data[left_child]) <= 0)
                {
                    break;
                }

                //Swap parent and child
                T tmp = data[parent_idx];
                data[parent_idx] = data[left_child];
                data[left_child] = tmp;
                parent_idx = left_child;

            }

            return front_item;
        }

        public int Count()
        {
            return data.Count;
        }

    }

}
