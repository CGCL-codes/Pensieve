#ifndef __MEM_H
#define __MEM_H
#include <stdio.h>  
#include <unistd.h>  
#include <sys/time.h>  
#include <string.h>  
#include <stdlib.h>  
#include <error.h>
  
#define VMRSS_LINE 16

//获取进程占用内存  
unsigned int get_proc_mem(unsigned int pid){  
      
    char file_name[64]={0};  
    FILE *fd;  
    char line_buff[512]={0};  
    sprintf(file_name,"/proc/%d/status",pid);
      
    fd =fopen(file_name,"r");  
    if(nullptr == fd){  
        return 0;  
    }  
      
    char name[64];  
    int vmrss;  
    for (int i=0; i<VMRSS_LINE-1;i++){  
        fgets(line_buff,sizeof(line_buff),fd);  
    }  

    fgets(line_buff,sizeof(line_buff),fd);
    sscanf(line_buff,"%s %d",name, &vmrss);
    fclose(fd);  
    // printf("%s\n", line_buff);  
    return vmrss;  
}  
 
//进程本身  
int get_pid(const char* process_name)  
{  
    char cmd[512];  
    sprintf(cmd, "pgrep %s -u tangwei", process_name);   
    
    FILE *pstr = popen(cmd, "r");
    
    if(pstr == nullptr){  
        printf("cmd error code : %d\n", errno);
        return 0;
    }

    char buff[512];  
    ::memset(buff, 0, sizeof(buff));  
    if(NULL == fgets(buff, 512, pstr)){  
        return 0;  
    }  
  
    return atoi(buff);  
} 


void print_mem(const char * process_name) {
    // printf("%s\n", process_name);
    int pid = get_pid(process_name);
    // printf("pid:%d\n", pid);
    malloc_trim(0);
    printf("current mem using : %u kb\n", get_proc_mem(pid));
}

void print_mem_without_trim(const char * process_name) {
    int pid = get_pid(process_name);
    printf("current mem using : %u kb\n", get_proc_mem(pid));
}

void print_mem(const int pid) {
    malloc_trim(0);
    printf("current mem using : %u kb\n", get_proc_mem(pid));
}

void print_mem_without_trim(const int pid) {
    printf("current mem using : %u kb\n", get_proc_mem(pid));
}

#endif
