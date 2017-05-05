		    # include<stdio.h>
			# include<time.h>
			# include<stdlib.h>
			# include<math.h>
			
			typedef struct nodedata{
				            int val;
			                struct nodedata *next;
			               }node;  
			
			
			
			node* makeq (node*,int);
			node* deleteq (node*);
			void showq (node*);
			int poisson_generation(double);
			double diffusion_transport(int,int);
			int active_transport(int,int);
			int subgrad_receiver(int [], int [], int,int); 
			double min(double,double);
			int transmitter(int);
			int metric_calculation();
			
			node* start[25];
			node *temp;
			int	gap_index,i,arr_pkt,sw_car,sw_pkt,src,des,src_n,b[25],count,count2;
			int tot_arv_calcium=0,tot_arv_BDNF=0,tot_lost=0,tot_dep=0,tot_wait=0,value,Level;
			double mean_arr_rate=20,mean_srv_rate= 60,l_mt=8*pow(10,-9),v_cargo=800*pow(10,-9),mean_sw_rate;
			double mean_arr_rate_ca=40,mean_srv_rate_ca= 80;
			double pr[3];
			int wt[]={31,38,25,30,48,26};	  // Quantity of each gate	
		    int val[]={25,21,18,29,35,23}; 	 // time taken
		    int W;							 // total packet served
		    int n=6; 						//No. of cargo
			int number_of_cargoes,m_to_m[3];
			int x;	
			int loss_subgrad;
			double fk;
			double x_distance=pow(10,-9),y_distance=0,Diffusion_coefficient=5.3*pow(10,-10),time_interval=pow(10,-9);
			
			
			
			int active_transport(int num_calcium, int num_BDNF)
			{
			// source to first level of tubules ( fick's law)     
			   	src_n =0 ;
			      for(src_n=0;src_n<3;src_n++)
			      {
				    tot_arv_calcium=tot_arv_calcium+num_calcium;
				    tot_arv_BDNF = tot_arv_BDNF+num_BDNF;
			        for(arr_pkt=1;arr_pkt<=num_BDNF;arr_pkt++)
			           {
			              count=0;
			              temp= start[src_n];     // checking ip q status
			              while(temp!=NULL)
			               {
			                 count++;
			                 temp=temp->next;
			               }
			           
			                if(count<b[src_n])    // if buffer is not full 
			               {
			               	  value=1;
			                  start[src_n]= makeq(start[src_n],value); // entry in ip q 
			                }
			                  if(count>=b[src_n])           // if buffer is full
			                    tot_lost++;  // packet drop
			            }      
				   }
				  
				    // cargoes taking ions on microtubule and storing it on other side queues
				    double calcium_concentration=0;
			    for(Level=0;Level<4;Level++)
				  {
				    des=6*Level+2;
				  for(src=6*Level;src<6*Level+3;src++)
				   {
				     des++;
				     
				     calcium_concentration = diffusion_transport(src*10,tot_arv_calcium);
				     mean_sw_rate= (v_cargo/l_mt)-10*calcium_concentration; 
				     number_of_cargoes= mean_sw_rate*0.9+(mean_sw_rate*0.1)*(rand()%1000)/1000;
			         for(sw_car=1;sw_car<=number_of_cargoes;sw_car++)
			         {
			         	for(sw_pkt=1;sw_pkt<=70;sw_pkt++)
			         	{
			         	    count=0;
			                temp= start[src];
			                
			               if(start[src]==NULL)
			                 break;
			               if(start[src]!=NULL)                // ip q not empty
			               {
			               	
			            	start[src] = deleteq(start[src]);
			            	temp= start[des];            // checking op q status
			            	
			                while(temp!=NULL)
			                {
			                 	
			                    count++;
			                    temp=temp->next;
			                }
			            	
			                if(count<b[des])     // if buffer not full
			                  {
			                  	value=1;
			                     start[des]= makeq(start[des],value);  // entry in op q         
			                  }
			                if(count>=b[des])            // if buffer full
			                     tot_lost++;   // packet drop
			            		
					     	}	
					    } 
			          }
			       }  
			       
			       if(Level==3)
			         continue;
			       
			       // tubule to tubule (to be added proper mathematical expression, now just using some probability) 
			       pr[0]=0.5,pr[1]=0.3,pr[2]=0.2; 				   
			       
			      for(src=6*Level+3;src<6*Level+6;src++)
			       { 
			         count = 0;
			         temp=start[src];
			         while(temp!=NULL)
			         {
			         	count++;
			         	temp=temp->next;
					 }
			         
			            for(des=6*Level+6;des<6*Level+9;des++)
			            {
			               	gap_index=(des-src)%3;
			            	m_to_m[gap_index]= (int)(count * pr[gap_index]);
			            	for(sw_pkt=1;sw_pkt<=m_to_m[gap_index];sw_pkt++)
			         	     {
			         	  	    // remove packet from that
			            	    start[src] = deleteq(start[src]);
			            	
			            		count2=0;
			            	 	temp= start[des];
			               
			                 	   while(temp!=NULL)
			                      {
			                        count2++;
			                        temp=temp->next;
			                      }
			            
			                 if(count2<b[des])     // if buffer not full
			                  {
			                  	value=1;
			                     start[des]= makeq(start[des],value);  // entry in op q           
			                  }
			                if(count2>=b[des])            // if buffer full
			                     tot_lost++;   // packet drop
			              
						   }
			
						}
			             
				   }
				   
		       for(i=des-2;i<=des;i++)
			         {
			         	while(start[i]!=NULL)
			         	{
			         	   start[i] = deleteq(start[i]);
						   tot_dep++;	
						 }
			         	     
					 }
			   	
			    }
			         
			   	    
			      
				    for(i=0;i<=des-3;i++)
			         {
			         	temp=start[i];
			         	while(temp!=NULL)
			         	 {
			         	 	tot_wait++;
			         	 	temp=temp->next;
						  }
					 }  
				 
			}
			
			double diffusion_transport(int distance,int initial_conc)
			{
				double D = 1 * pow(10,-9);
				double B = 0.5;
				double x = distance*pow(10,-9);
				double time =  1000* pow(10,-9);
				double Deff = D/(1+B);
				double expo_term = -x*x/(4*Deff*time);
				double constant_term = 1/(pow((4*3.142*Deff*time),1.5));
				double value = initial_conc * constant_term*exp(expo_term)*(pow(10,-27));
				return value;
			}
				
				
				int poisson_generation(double lamda)
			{
				double L,p=1;
				int k=0;
				double u;
				L=1/exp(lamda);
				do
				{
					k=k+1;
					u=rand()%10000;
					u=u/10000;
					p=p*u;
				} while(p>L);
			
				return(k-1);
				
			}
			
			
			int transmitter(int arrival_rate)
			{
				int number_of_elements_generated = poisson_generation(arrival_rate)*pow(10,-6)*sqrt(pow(x_distance,2)+pow(y_distance,2))/(2*sqrt(Diffusion_coefficient*3.14)*pow(time_interval,1.5));;
				return number_of_elements_generated;	
			}
						
			int main()
			{
			   time_t t;
			   srand((unsigned)time(&t));
			   
			   for(int loop=0;loop<25;loop++)
			   {
			    start[loop]=NULL;
			    b[loop]=11000;
			   }
			   for(int T=0;T<15;T++)
			   {
			   	    int number_of_calcium_generated= transmitter(mean_arr_rate_ca);
					int number_of_BDNF_generated = transmitter(mean_arr_rate);
			   		int trans = active_transport(number_of_calcium_generated, number_of_BDNF_generated);
			   }
			   		W= tot_dep;
					printf("\nTime taken ::: %d unit",subgrad_receiver(val,wt,W,n));
			   		int metrics = metric_calculation();
			   		return 0;
				
			}
			
			int metric_calculation()
			{
				double dia_MT = 4*pow(10,-9);
			    double area_of_cross_section_of_MT = 3.142*pow(dia_MT,2)/4;
			    double rate_MT = 1*pow(10,-6)*70*pow(10,-6)*(pow(10,-9))/v_cargo;
			    double flux_MT = rate_MT/area_of_cross_section_of_MT;
			  	
		     	double c0 =poisson_generation(mean_arr_rate)*pow(10,-6);
				double r11= sqrt((x_distance*x_distance)+(y_distance*y_distance));
				double Jx = c0*r11*exp(-r11*r11/(4*Diffusion_coefficient*time_interval))/(4*time_interval*sqrt(Diffusion_coefficient*time_interval*3.142));
				double Jy = c0*r11*exp(-y_distance*y_distance/(4*Diffusion_coefficient*time_interval))/(4*time_interval*sqrt(Diffusion_coefficient*time_interval*3.142));
				double delay_diffusion = (x_distance*x_distance + y_distance*y_distance)/(2.0*Diffusion_coefficient);
				double delay_MT_track = l_mt/v_cargo * 4;
				 
			   
			   printf("\ntotal BDNF molecules arrived: %d",tot_arv_BDNF);
			   printf("\ntotal BDNF molecules served: %d",tot_dep);
			   printf("\ntotal BDNF molecules waiting in inner stage queues now: %d",tot_wait);
			   printf("\n\noptimised number of gates for the served BDNF molecules.\n");
			   
			   
			   printf("\ntotal BDNF molecules lost in MT: %d\n",tot_lost);
			   printf("\ntotal BDNF molecules lost in optimisation part: %d\n",loss_subgrad);
			   printf("\ntotal BDNF molecules lost in whole model: %d\n",tot_lost+loss_subgrad);
			   printf("\nblocking probability= %lf percent \n",100*((tot_lost+loss_subgrad)*1.0/tot_arv_BDNF));
			   printf("\nMessage Deliverability = %lf percent \n",100-(100*((tot_lost+loss_subgrad)*1.0/tot_arv_BDNF)));
			   printf("\nfunction = %g \n", fk );
               printf("\naverage delay in clearance = %g \n",1/(fk*W)); 
               printf("\nDelay due to diffusion = %g\n",delay_diffusion);
               printf("\nDelay due on MT tracks = %g\n",delay_MT_track);
               printf("\nDiffusive flux in x direction = %g \n",Jx);
			   printf("\nDiffusive flux in y direction = %g \n",Jy);
			   printf("\nDiffusive flux = %g \n",Jx+Jy);
			   printf("\nFlux on MT = %g  area = %g  rate = %g\n",flux_MT,area_of_cross_section_of_MT,rate_MT);
			   return 0;
			}
			
				node* makeq (node* start1,int value)
			{
			    node* t=(node*)malloc(sizeof(node));
			    node* t1;
			   t->val= value;
			   
			      if(start1==NULL)
			         {
			           start1=t;
			           t1=t;
			           t1->next=NULL;
			           return (start1);
			         }
			      else
			         {
			            t1=start1;
			            
			               while(t1->next!=NULL)
			                {
			                  t1=t1->next;
			                }
			
			           t1->next=t;
			           
			           t->next=NULL;
			           return(start1);
			        
			         } 
			         
			    
			}
			
			node* deleteq (node* start1)
			{
			   node* t,*t1, *last;
			   int a,b;
			  
			
			     t1=start1;
			          while(t1->next!=NULL)
			                t1=t1->next;
		
			     if(start1==t1)
			      {
			       t=start1;
			       a=t->val;
			       x=a;
			       free(t);
			       start1=NULL;
			       
			       //printf("\ndeleted value=%d",a);
			       return(start1);
			      }
			    
			     else
			      {
			       t=start1;
			       a=t->val;
			      
			       x=a;
			       start1=t->next;
			       free(t);
			      // printf("\ndeleted value=%d",b);
			       return(start1);
			    
			      }
			     
			}
			
			void showq (node* start1)
			{
			  node* t=start1;
			    if(start1==NULL)
			     printf("\nqueue is empty");
			    else
			   {
			    while(t!=NULL)
			     {
			        printf("%d\n",t->val);
			     
			        t=t->next;
			     }
			    }
			}
			

 int subgrad_receiver(int val[], int wt[], int W,int n)
{
	double zip=pow(10,20);
	double lambda1 = 20;
	double lambda2 = 120;
	double lambda3 = 14;
	double lambda4 = 33;
	double lambda5 = 1122;
	double lambda6 = 125;
	int k=0,imp=0;
	double best_lambda1=lambda1;
	double best_lambda2=lambda2;
	double best_lambda3=lambda3;
    double best_lambda4=lambda4;
    double best_lambda5=lambda5;
	double best_lambda6=lambda6;
	double best_L_lambda=-10000000;
	float delta = 0.7;
	double x1,x2,x3,x4,x5,x6;
	double c1,c2,c3,c4,c5,c6;
	int a;
	double r1= val[0],r2=val[1],r3=val[2],r4=val[3],r5=val[4],r6=val[5];
	c1=wt[0],c2=wt[1],c3=wt[2],c4=wt[3],c5=wt[4],c6=wt[5];
	double L_lambda=0,fk,tk,G1,G2,G3,G4,G5,G6;
	tk=10;
	a=2;
	
	do
   {
	abcd :
		
	imp++;
	k++;
	

	x1=pow(((r1*lambda1*(c1+c2+c3+c4+c5+c6)/(W*c1)))*W/(r1*a),1/(a-1))*W/r1;
	x2=pow(((r2*lambda2*(c1+c2+c3+c4+c5+c6)/(W*c2)))*W/(r2*a),1/(a-1))*W/r2;
	x3=pow(((r3*lambda3*(c1+c2+c3+c4+c5+c6)/(W*c3)))*W/(r3*a),1/(a-1))*W/r3;
	x4=pow(((r4*lambda4*(c1+c2+c3+c4+c5+c6)/(W*c4)))*W/(r4*a),1/(a-1))*W/r4;
	x5=pow(((r5*lambda5*(c1+c2+c3+c4+c5+c6)/(W*c5)))*W/(r5*a),1/(a-1))*W/r5;
	x6=pow(((r6*lambda6*(c1+c2+c3+c4+c5+c6)/(W*c6)))*W/(r6*a),1/(a-1))*W/r6;
	
	L_lambda=pow((r1*x1)/W,a)+pow((r2*x2)/W,a)+pow((r3*x3)/W,a)+pow((r4*x4)/W,a)+pow((r5*x5)/W,a)+pow((r6*x6)/W,a)+lambda1*(1-(r1*x1*(c1+c2+c3+c4+c5+c6)/(W*c1)))+lambda2*(1-(r2*x2*(c1+c2+c3+c4+c5+c6)/(W*c2)))+lambda3*(1-(r3*x3*(c1+c2+c3+c4+c5+c6)/(W*c3)))+lambda4*(1-(r4*x4*(c1+c2+c3+c4+c5+c6)/(W*c4)))+lambda5*(1-(r5*x5*(c1+c2+c3+c4+c5+c6)/(W*c5)))+lambda6*(1-(r6*x6*(c1+c2+c3+c4+c5+c6)/(W*c6)));//+lambda3*((-c1*x1/W-c2*x2/W)+1);
	fk=pow((r1*x1)/W,a)+pow((r2*x2)/W,a)+pow((r3*x3)/W,a)+pow((r4*x4)/W,a)+pow((r5*x5)/W,a)+pow((r6*x6)/W,a);
	
	if(L_lambda>best_L_lambda)
	{
		best_L_lambda=L_lambda;
		best_lambda1=lambda1;
		best_lambda2=lambda2;
		best_lambda3=lambda3;
		best_lambda4=lambda4;
		best_lambda5=lambda5;
		best_lambda6=lambda6;
		imp=-1;
	}
	

	if(x1>0 && x2>0&& x3>0&& x4>0 && x5>0 &&x6>0)
	{
		if(fk<zip)
		zip=fk;
    }
   

	if(imp>20)
	{
		delta=delta/2;
		lambda1=best_lambda1;
		lambda2=best_lambda2;
		lambda3=best_lambda3;
		lambda4=best_lambda4;
		lambda5=best_lambda5;
		lambda6=best_lambda6;
		imp=0;
		goto abcd;
	}
	

	double error= abs((zip-best_L_lambda)/best_L_lambda);
	if(k>1000000 || delta < 0.0000001 || tk< 0.00000001 || error < 0.00000001)
	   break;
	
	G6 = 1-(r6*x6*(c1+c2+c3+c4+c5+c6)/(W*c6));
	G5 = 1-(r5*x5*(c1+c2+c3+c4+c5+c6)/(W*c5));
	G4 = 1-(r4*x4*(c1+c2+c3+c4+c5+c6)/(W*c4));
	G3 = 1-(r3*x3*(c1+c2+c3+c4+c5+c6)/(W*c3));
	G2 = 1-(r2*x2*(c1+c2+c3+c4+c5+c6)/(W*c2));
	G1 = 1-(r1*x1*(c1+c2+c3+c4+c5+c6)/(W*c1));
	tk= delta * (zip-L_lambda)/(G1*G1+G2*G2+G3*G3+G4*G4+G5*G5+G6*G6);
    
    
	lambda1=min(-1,lambda1+(tk*G1));
    lambda2=min(-1,lambda2+(tk*G2));
    lambda3=min(-1,lambda3+(tk*G3));
    lambda4=min(-1,lambda4+(tk*G4));
    lambda5=min(-1,lambda5+(tk*G5));
	lambda6=min(-1,lambda6+(tk*G6));
	}while(true);


printf("\n No. of gates of capacity %f required : %f \n",c1,ceil(x1));
printf("\n No. of gates of capacity %f required : %f \n",c2,ceil(x2));
printf("\n No. of gates of capacity %f required : %f \n",c3,ceil(x3));
printf("\n No. of gates of capacity %f required : %f \n",c4,ceil(x4));
printf("\n No. of gates of capacity %f required : %f \n",c5,ceil(x5));
printf("\n No. of gates of capacity %f required : %f \n",c6,ceil(x6));

loss_subgrad = W-(c1*ceil(x1)+c2*ceil(x2)+c3*ceil(x3)+c4*ceil(x4)+c5*ceil(x5)+c6*ceil(x6));

fk=pow((r1*x1)/W,a)+pow((r2*x2)/W,a)+pow((r3*x3)/W,a)+pow((r4*x4)/W,a)+pow((r5*x5)/W,a)+pow((r6*x6)/W,a);

printf("%g\n",fk);

return 0;
    }
    
    
	    double min(double a,double b)
	   {
		if(a<b)
		return b;
		else
		return a;
		}
