#include <bits/stdc++.h>
using namespace std;
double eyeX, eyeY, eyeZ, lookX, lookY, lookZ, upX, upY, upZ, fovY, aspectRatio, near, far;
stack<int> pushed_idx;
stack < vector< vector<double>> > commands_stack;
string stage1 ,stage2,stage3;
void readFirst4Lines(ifstream *myfile, string line)
{
    //***********  line1  ************
    getline(*myfile, line);
    stringstream stream(line);
    string temp;

    for(int i=0; i<3; i++)
    {
        getline(stream,temp,' ');
        if(i==0) eyeX = stod(temp);
        else if(i==1) eyeY = stod(temp);
        else if(i==2) eyeZ = stod(temp);
    }
    cout << eyeX << " "<< eyeY << " " << eyeZ << endl;

    // line 2
    getline(*myfile, line);
    stringstream stream1(line);
    //string temp;

    for(int i=0; i<3; i++)
    {
        getline(stream1,temp,' ');
        if(i==0) lookX = stod(temp);
        else if(i==1) lookY= stod(temp);
        else if(i==2) lookZ = stod(temp);
    }
    cout << lookX << " " << lookY << " " << lookZ << endl;


    //    line 3
    getline(*myfile,line);
    stringstream stream2(line);

    for(int i=0; i<3; i++)
    {
        getline(stream2,temp,' ');
        if(i==0) upX = stod(temp);
        else if(i==1) upY= stod(temp);
        else if(i==2) upZ= stod(temp);

    }
    cout << upX << " " << upY<< " " << upZ<< endl;


// line  4
    getline(*myfile,line);
    stringstream stream3(line);

    for(int i=0; i<4; i++)
    {
        getline(stream3,temp,' ');
        if(i==0) fovY = stod(temp);
        else if(i==1) aspectRatio = stod(temp);
        else if(i==2) near = stod(temp);
        else if(i==3)   far = stod(temp);
    }

}

vector<double> add_vectors(vector<double> a, vector<double>b)
{
    vector<double> result(3);
    result[0] = a[0] + b[0];
    result[1] = a[1] + b[1];
    result[2] = a[2] + b[2];
    return result;
}

vector<double> scalerMult_vector(double n, vector<double>b)
{
    vector<double> result(3);
    result[0] = round(b[0]*n);
    result[1] = round(b[1]*n);
    result[2] = round(b[2]*n);
    return result;
}

double dot_product(vector<double> a, vector<double>b)
{
    double result;
    result += a[0] * b[0];
    result += a[1] * b[1];
    result += a[2] * b[2];
    return result;
}


vector<double> cross_product(vector<double> a, vector<double>b)
{
    vector<double> result(3);
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];

    return result;
}

vector <vector<double>>  multiplyVector( vector< vector<double>> Mat1, vector< vector<double>> Mat2 )
{
    int n = Mat1.size();     // a rows
    int m = Mat1[0].size();  // a cols
    int p = Mat2[0].size();  // b cols

    vector <vector<double>> result(n, vector<double>(p, 0));

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < p; ++j)
        {
            for (int k = 0; k < m; ++k)
            {
                result[i][j] += Mat1[i][k] * Mat2[k][j];
            }
        }
    }
    return result;

}


vector<double> R_formula(vector<double> axis_unit_vect , vector<double> a , double angle)
{
    vector<double> vect1 = scalerMult_vector( cos(angle*3.1416/180) ,axis_unit_vect);
    vector<double> vect2 = scalerMult_vector((1-cos(angle*3.1416/180)) * dot_product(a,axis_unit_vect), a   );
    vector<double> vect3 = scalerMult_vector(sin((angle*3.1416/180)) , cross_product(a,axis_unit_vect)  );

    vector<double> result = add_vectors(vect1,vect2);
    result = add_vectors(result,vect3);
    return result;
}


void print4_4Vector(vector<vector<double>> arg)
{
    cout << "Matrix: " << endl;
    for(int i=0; i<arg.size(); i++)
    {
        for(int j=0; j<arg[0].size(); j++)
        {
            cout << arg[i][j] << " ";

        }
        cout << endl;
    }
}



void translate_command(double tx,double ty,double tz)
{
    vector<vector<double> >  transtaror_matrix(4,vector<double>(4,0)) ;
    for(int i=0; i<4; i++)
    {
        transtaror_matrix[i][i] = 1;
    }
    transtaror_matrix[0][3] = tx;
    transtaror_matrix[1][3] = ty;
    transtaror_matrix[2][3] = tz;
    //cout << "translator\n" ;
    //print4_4Vector(transtaror_matrix);
   // cout << "stacktop ";
    //print4_4Vector(commands_stack.top());
    commands_stack.push(multiplyVector(commands_stack.top(),transtaror_matrix));
    //cout << "stacktop ";
   // print4_4Vector(commands_stack.top());
}






void rotate_command( double angle , double ax , double ay , double az)
{
    double normalizer = sqrt(ax*ax + ay*ay + az*az);
    vector<double> a ;
    a.push_back(ax/normalizer);
    a.push_back(ay/normalizer);
    a.push_back(az/normalizer);

    vector<double> i = {{1.0}, {0.0}, {0.0} };
    vector<double> j = {{0.0}, {1.0}, {0.0} };
    vector<double> k = {{0.0}, {0.0}, {1.0} };


    vector<double> c1 = R_formula(i,a,angle);
    vector<double> c2 = R_formula(j,a,angle);
    vector<double> c3 = R_formula(k,a,angle);

    vector<vector<double>> rotation_matrix = { { { c1[0] }, { c2[0] }, { c3[0] } , {0} },
                                                 { { c1[1] }, { c2[1] }, { c3[1] } , {0} },
                                                 { { c1[2] }, { c2[2] }, { c3[2] } , {0} },
                                                 {  {0},        {0},        {0},    {1.0} }
                                                };


    //cout << "rotation" ;
    //print4_4Vector(rotation_matrix);
    commands_stack.push(  multiplyVector(commands_stack.top(),rotation_matrix));
   // cout << "rotationtop";
    //print4_4Vector(commands_stack.top());
}

void scale_command(double sx, double sy,double sz)
{
     vector<vector<double> >  scale_matrix(4,vector<double>(4,0)) ;
     scale_matrix[0][0] =sx;
     scale_matrix[1][1] = sy;
     scale_matrix[2][2] = sz;
     scale_matrix[3][3] = 1;
     commands_stack.push(multiplyVector(commands_stack.top(),scale_matrix));
     //cout << "Scale:";
     //print4_4Vector(commands_stack.top());
}



vector<vector<double>> scalePoint( vector<vector<double>> point ){
    if(point[3][0] != 1 ){
        for( int i = 0 ; i < 4; i++ ){
            point[i][0] = point[i][0] / point[3][0];
        }
    }

    return point;
}


vector<double>  normalize(vector<double> c)
{
     double val = 0.0;

     for(int i=0; i<c.size() ; i++) {
        val += c[i]*c[i];
     }
     val = sqrt(val);
     c[0] = c[0]/val;
     c[1] = c[1] /val;
     c[2] = c[2] / val;
     return c;
}

vector<vector<double> > viewTransform()
{
    vector<vector <double>> T_vect = { {  {1.0}, {0.0},  {0.0} ,  {-eyeX}  },
                                        {  {0.0}, {1.0},  {0.0} ,   {-eyeY}  } ,
                                        {  {0.0}, {0.0},  {1.0} ,   {-eyeZ}  } ,
                                        {  {0.0}, {0.0},  {0.0} ,   {1.0}    }
                                        };


    vector<double> L = { {lookX - eyeX}, {lookY - eyeY} , {lookZ- eyeZ} } ;
    L = normalize(L);
    vector<double> up;
    up.push_back(upX);
    up.push_back(upY);
    up.push_back(upZ);
    vector<double> R = cross_product(L,up);
    R = normalize(R);
    vector< double > U = cross_product(R,L);
    U = normalize(U);
    vector<vector <double>> R_vect = {
                                        { {R[0]},   {R[1]},     { R[2] } ,  { 0.0 }    },
                                        { {U[0]},   {U[1]},     { U[2] } ,  { 0.0 }    },
                                        { {-L[0]},  {-L[1]},    {- L[2]} , { 0.0 }    },
                                        { { 0.0 },  { 0.0 },    { 0.0 } ,   { 1.0 }    },

                                        } ;
     cout << "Lvect " << L[0]  << " " << L[1] << " " << L[2] << endl;

    return multiplyVector(R_vect, T_vect);


}



vector<vector<double>> projectionTransform()
{
    double fovX = fovY * aspectRatio;
    double t = near * tan((fovY/2) * 3.1416/180);
    double r = near * tan((fovX/2) * 3.1416/180);

    return {
            { {near/r}, {0.0} , {0.0} , {0.0} },
            { {0.0}, {near/t},   {0.0} , {0.0} },
            { {0.0} , {0.0}  ,  { -(far+near)/(far-near) },  {   -(2*far*near)/(far-near) }  },
            {  {0.0} , {0.0} , {-1.0} , {0.0}     }

        }  ;

}



void TransformPoint(vector<vector<double>>  point )
{
    point = multiplyVector(commands_stack.top(), point);
    //cout << "transformed point= ";
    //print4_4Vector(point);
    point = scalePoint(point);

    for(int i=0; i<3; i++)
    {
        stage1 +=to_string( point[i][0] ) + "  ";
    }
    stage1 += "\n";


    point = multiplyVector(viewTransform() , point);
    cout << "print disi\n";
    print4_4Vector(point);
    point = scalePoint(point);
    print4_4Vector(point);

    for(int i=0; i<3; i++)
    {
        stage2 +=to_string( point[i][0] ) + "  ";
    }
    stage2 += "\n";


    point = multiplyVector(projectionTransform(),point);
    point = scalePoint(point);

    for(int i=0; i<3; i++)
    {
        stage3 +=to_string( point[i][0] ) + "  ";
    }
    stage3 += "\n";


}

vector <vector<double>> readPoint(ifstream *myfile, string line)
{
    getline(*myfile,line);
    stringstream stream(line);
    string temp;
    vector<vector<double>> point(4,vector<double>(1,0));
    for(int i=0; i<3; i++)
    {
        getline(stream,temp,' ');
        cout << temp << endl;
        point[i][0] = stod(temp);
    }

    point[3][0] = 1.0;
    return point;
}

int main()
{
    ifstream myfile;
    string line,command;
    myfile.open("scene.txt");

    readFirst4Lines(&myfile,line);

    vector <vector<double>> firstElement(4,vector<double>(4,0));
    for(int i=0; i<4; i++)
    {
        firstElement[i][i] = 1.0;
    }


    commands_stack.push(firstElement);


    print4_4Vector(commands_stack.top());


    while (getline(myfile, command))
    {
        cout << command << endl;

        if(command == "triangle")
        {
            cout << "tri paise\n";
            for(int i=0; i<3; i++)
            {
                vector<vector<double>> point1 = readPoint(&myfile,line);
                TransformPoint(point1);
                print4_4Vector(point1);
            }

            stage1 += "\n";
            stage2 += "\n";
            stage3 += "\n";
        }

        else if(command == "translate")
        {
            cout << "translate paise\n";
            getline(myfile,line);
            stringstream stream(line);
            string temp;
            double tx,ty,tz;
            for(int i=0; i<3; i++)
            {
                getline(stream,temp,' ');
                cout << temp << endl;
                if(i==0) tx =  stod(temp);
                else if(i==1) ty =  stod(temp);
                else if(i==2) tz =  stod(temp);
            }

            translate_command(tx,ty,tz);
            print4_4Vector(commands_stack.top());

        }

        else  if(command == "scale")
        {
            cout << "scale paise\n";
            getline(myfile,line);
            stringstream stream(line);
            string temp;
            double sx,sy,sz;
            for(int i=0; i<3; i++)
            {
                getline(stream,temp,' ');
                cout << temp << endl;
                if(i==0) sx =  stod(temp);
                else if(i==1) sy =  stod(temp);
                else if(i==2) sz =  stod(temp);
            }

            scale_command(sx,sy,sz);
            print4_4Vector(commands_stack.top());
        }

        else  if(command == "rotate")
        {
            cout << " Rotate pse\n";

            getline(myfile,line);
            stringstream stream(line);
            string temp;
            double angle , ax,ay,az;
            for(int i=0; i<4; i++)
            {
                getline(stream,temp,' ');
                cout << temp << endl;
                if(i==0) angle = stod(temp);
                else if(i==1) ax =  stod(temp);
                else if(i==2) ay =  stod(temp);
                else if(i==3) az =  stod(temp);
            }

            rotate_command(angle,ax,ay,az);
            print4_4Vector(commands_stack.top());

        }

        else  if(command == "push")
        {
            cout << "pusssss\n";
            pushed_idx.push(commands_stack.size()-1);
            print4_4Vector(commands_stack.top());
        }
        else  if(command == "pop")
        {
            cout << "pooooop\n";
            if(commands_stack.size() > 1)
            {
                int idx = pushed_idx.top();
                pushed_idx.pop();
                int noOfmatricesTobePoped = commands_stack.size()- idx - 1;

                for(int i=0; i<noOfmatricesTobePoped; i++)
                {
                    commands_stack.pop();
                }

            }
            print4_4Vector(commands_stack.top());

        }

        else  if(command == "end")
        {
            cout << "sesh\n";
            break;
        }
    }
    cout << "kaj korse";
    cout << "stage1 " << stage1;
    cout << endl << "stage 2 \n";
    cout << stage2;
    cout << endl << "stage 3 \n";
    cout << stage3;
    ofstream teststage1,teststage2,teststage3;
    teststage1.open("testStage1.txt");
    teststage2.open("testStage2.txt");
    teststage3.open("testStage3.txt");
    teststage1 << stage1;
    teststage2 << stage2;
    teststage3 << stage3;
    teststage1.close();
    teststage2.close();
    teststage3.close();
    return 0;
}
