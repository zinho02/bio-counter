#include <stdio.h>

int main(){
	int ini;
	double soma = 0.1;
	while(scanf("%d", &ini)!= EOF){
		
		printf("if (i > %d && i < %d)\n",ini, ini+600 );
		printf("  pulso(i) = pulso(i-1)+%lf;\n",soma);
		printf("endif\n");
		
		printf("if (i >= %d && i < %d)\n", ini+600,ini+900);
		printf("  pulso(i) = pulso(i-1)-%lf;\n",soma*2);
		printf("endif\n");
	}
	return 0;
}